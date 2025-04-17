const dotenv = require('dotenv');
dotenv.config();



const express = require("express");
const { Pool } = require('pg');
const { v4: uuidv4 } = require("uuid");
const fs = require("fs");


const app = express();
const port = process.env.PORT || 3000;
const MAX_PROTEIN_LENGTH = process.env.MAX_PROTEIN_LENGTH || 2000; 
//const DATA_PATH = "./data/proteins.json";
app.use('/files', express.static('data'));
app.use(express.json());
app.use(express.urlencoded({ extended: true }));
app.use(express.text());
const WINDOW_SIZE = 17;
// Create connection pool
const pool = new Pool({
    host: process.env.PG_HOST,
    port: process.env.PG_PORT,
    database: process.env.PG_DATABASE,
    user: process.env.PG_USER,
    password: process.env.PG_PASSWORD
    });

    
    // Test the database connection
    pool.query('SELECT NOW()', (err, res) => {
        if (err) {
            console.error('Error connecting to PostgreSQL:', err);
            process.exit(1);
        } else {
            console.log('Connected to PostgreSQL:', res.rows[0]);
            app.listen(port, () => {
                console.log(` Server running on port ${port}`);
            });
        }
    });
    

let proteinData = null; 

const VALID_AMINO_ACIDS = /^[ACDEFGHIKLMNPQRSTVWY]+$/;
const AMINO_ACID_WEIGHTS = {
    A: 89.09, R: 174.2, N: 132.12, D: 133.1, C: 121.16,
    E: 147.13, Q: 146.15, G: 75.07, H: 155.16, I: 131.18,
    L: 131.18, K: 146.19, M: 149.21, F: 165.19, P: 115.13,
    S: 105.09, T: 119.12, W: 204.23, Y: 181.19, V: 117.15
};

function generateProteinId() {
    return uuidv4();
}

function readData(filePath) {
    try {
        return fs.existsSync(filePath) ? JSON.parse(fs.readFileSync(filePath)) : { proteins: [] };
    } catch (error) {
        console.error(`Error reading file ${filePath}:`, error);
        return { proteins: [] }; 
    }
}


function writeData(filePath, data) {
    fs.writeFileSync(filePath, JSON.stringify(data, null, 2));
}

function calculateMolecularWeight(sequence) {
    return sequence.split("")
        .reduce((total, aa) => total + (AMINO_ACID_WEIGHTS[aa] || 0), 0)
        .toFixed(2);
}
function populateExtraFields(proteinData) {
    if (!proteinData.molecularWeight && proteinData.sequence) {
        proteinData.molecularWeight = calculateMolecularWeight(proteinData.sequence);
    }
    return proteinData;
}

function getCurrentTimestamp() {
    return new Date().toISOString().split(".")[0] + "Z";
}



class NotFoundError extends Error {
    constructor(message) {
        super(message);
        this.name = "NotFoundError";
    }
}

class ConflictError extends Error {
    constructor(message) {
        super(message);
        this.name = "ConflictError";
    }
}


//-------------Authentication---------------------------------//
async function authenticateUser(req, res, next) {
  const userId = req.header("X-User-ID");
  if (!userId) {
    return res.status(401).json({ error: "Missing X-User-ID header" });
  }
  try {
    const result = await pool.query("SELECT * FROM users WHERE id = $1", [userId]);
    if (result.rows.length === 0) {
      return res.status(401).json({ error: "Invalid user" });
    }
    req.user = result.rows[0];
    next();
  } catch (error) {
    console.error("Authentication error:", error);
    res.status(500).json({ error: "Internal Server Error" });
  }
}
app.use("/api", authenticateUser);



const GOR_PROPENSITY = {
  A: { H: 1.42, E: 0.83, C: 0.56 },
  R: { H: 0.98, E: 0.93, C: 0.89 },
  N: { H: 0.67, E: 0.89, C: 0.94 },
  D: { H: 1.01, E: 0.54, C: 1.46 },
  C: { H: 0.70, E: 1.19, C: 0.94 },
  E: { H: 1.51, E: 0.37, C: 1.02 },
  Q: { H: 1.11, E: 1.10, C: 0.98 },
  G: { H: 0.57, E: 0.75, C: 1.31 },
  H: { H: 1.00, E: 0.87, C: 0.95 },
  I: { H: 1.08, E: 1.60, C: 0.47 },
  L: { H: 1.21, E: 1.30, C: 0.59 },
  K: { H: 1.16, E: 0.74, C: 0.96 },
  M: { H: 1.45, E: 1.05, C: 0.60 },
  F: { H: 1.13, E: 1.38, C: 0.61 },
  P: { H: 0.57, E: 0.55, C: 1.52 },
  S: { H: 0.77, E: 0.75, C: 1.32 },
  T: { H: 0.83, E: 1.20, C: 0.96 },
  W: { H: 1.08, E: 1.37, C: 0.65 },
  Y: { H: 0.69, E: 1.47, C: 0.71 },
  V: { H: 1.06, E: 1.70, C: 0.48 },
};


function predictSecondaryStructure(proteinId, sequence) {
  let secondaryStructure = "";
  let confidenceScores = [];

  for (let i = 0; i < sequence.length; i++) {
      let likelihoods = { H: 1, E: 1, C: 1 };

      for (let j = -Math.floor(WINDOW_SIZE / 2); j <= Math.floor(WINDOW_SIZE / 2); j++) {
          let index = i + j;
          if (index >= 0 && index < sequence.length) {
              const aminoAcid = sequence[index];

              if (GOR_PROPENSITY[aminoAcid]) {
                  likelihoods.H *= GOR_PROPENSITY[aminoAcid]["H"];
                  likelihoods.E *= GOR_PROPENSITY[aminoAcid]["E"];
                  likelihoods.C *= GOR_PROPENSITY[aminoAcid]["C"];
              }
          }
      }
  
      let sorted = Object.entries(likelihoods).sort((a, b) => b[1] - a[1]);
      secondaryStructure += sorted[0][0]; 
      confidenceScores.push((sorted[0][1] - sorted[1][1]).toFixed(2)); 
  }

  return {proteinId, sequence, secondaryStructure, confidenceScores };
}


async function fragmentAndStoreSequence(client, proteinId, sequence) {
  if (
    !proteinId ||
    !/^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/.test(proteinId)
  ) {
    throw new Error("Invalid proteinId: Must be a valid UUID");
  }

  const windowSize = 15;
  const stepSize = 5;

  try {
    for (let i = 0; i <= sequence.length - windowSize; i += stepSize) {
      const fragment = sequence.substring(i, i + windowSize);
      const fragmentId = uuidv4();

    
      const { secondaryStructure, confidenceScores } = predictSecondaryStructure(proteinId, fragment);
      const motifs = identifyMotifs(fragment);

      // instead of JSON path
      const fragmentUrl = `http://localhost:3000/api/fragments/${fragmentId}`;

   
      await client.query(
        `INSERT INTO fragments (
          fragment_id, protein_id, sequence, start_position, end_position, secondary_structure, url
        ) VALUES ($1, $2, $3, $4, $5, $6, $7)`,
        [fragmentId, proteinId, fragment, i + 1, i + windowSize, secondaryStructure, fragmentUrl]
      );

      
      for (const motif of motifs) {
        await client.query(
          `INSERT INTO motifs (
            motif_id, fragment_id, motif_pattern, motif_type, start_position, end_position, confidence_score
          ) VALUES ($1, $2, $3, $4, $5, $6, $7)`,
          [
            uuidv4(),
            fragmentId,
            motif.pattern,
            motif.type,
            motif.start,
            motif.end,
            motif.confidence,
          ]
        );
      }
    }
  } catch (error) {
    console.error("Fragmentation error:", error);
    throw error;
  }
}




function identifyMotifs(fragment) {
  //console.log("Checking fragment for motifs:", fragment);  // Debugging log
  
  const motifs = [];
  const regexPatterns = [
    { type: "N-glycosylation site", pattern: /N[^P][ST][^P]/g },
    { type: "Casein kinase II phosphorylation site", pattern: /[ST].{2}[DE]/g },
    { type: "Tyrosine kinase phosphorylation site", pattern: /[RK].{0,2}[DE]/g }
  ];

  for (const motif of regexPatterns) {
    let match;
    while ((match = motif.pattern.exec(fragment)) !== null) {
    //console.log("Found Motif:", motif.type, "at index", match.index + 1, "in fragment:", fragment);

      motifs.push({
        type: motif.type,
        pattern: match[0],
        start: match.index + 1,
        end: match.index + match[0].length,
        confidence: Math.random().toFixed(2)
      });
    }
  }
  return motifs;
}

app.post("/api/proteins", async (req, res) => {
  const client = await pool.connect();

  try {
    const { sequence, name, description } = req.body;

   
    if (!sequence || sequence.length > 1000 || !/^[ACDEFGHIKLMNPQRSTVWY]+$/.test(sequence)) {
      return res.status(400).json({
        error: "Invalid sequence format or length exceeds 2000 characters."
      });
    }

    const generatedName = name && name.length <= 100
      ? name
      : `Protein_${sequence.slice(0, 8)}_${Math.floor(Date.now() / 1000)}`;

    const proteinId = uuidv4();
    const molecularWeight = calculateMolecularWeight(sequence);
    const sequenceLength = sequence.length;
    const createdAt = new Date().toISOString();
    const updatedAt = createdAt;


    const sequenceUrl = `http://localhost:3000/api/proteins/${proteinId}/sequence`;

    await client.query("BEGIN");

   
    await client.query(
      `INSERT INTO proteins (
        protein_id, name, description, molecular_weight, sequence_length, sequence_url, created_at, updated_at
      ) VALUES ($1, $2, $3, $4, $5, $6, $7, $8)`,
      [
        proteinId,
        generatedName,
        description || null,
        molecularWeight,
        sequenceLength,
        sequenceUrl,
        createdAt,
        updatedAt,
      ]
    );

    await fragmentAndStoreSequence(client, proteinId, sequence);

    await client.query("COMMIT");

    res.status(201).json({
      message: "Successfully created the protein",
      proteinId,
      name: generatedName,
      description,
      molecularWeight,
      sequenceLength,
      createdAt,
      updatedAt,
      sequenceUrl,
    });

  } catch (error) {
    await client.query("ROLLBACK");
    console.error("Transaction failed:", error);
    res.status(500).json({ error: "Internal Server Error" });
  } finally {
    client.release();
  }
});



//-----------------------------GET PROTIENS LIST--------------------------------------------------//
app.get("/api/proteins", async (req, res) => {
  const client = await pool.connect();
  try {
      
      let { limit = 10, offset = 0, sort } = req.query;

      const allowedParams = ['limit', 'offset', 'sort'];
      for (const param in req.query) {
        if (!allowedParams.includes(param)) {
          return res.status(400).json({ error: `Invalid query parameter: ${param}` });
        }
      }

      limit = parseInt(limit);
      offset = parseInt(offset);

     
      if (isNaN(limit) || isNaN(offset) || limit < 1 || offset < 0) {
          return res.status(400).json({ error: "Limit must be ≥ 1 and offset must be ≥ 0." });
      }

   
      const validSortFields = {
          name: "name",
          createdAt: "created_at",
          molecularWeight: "molecular_weight",
          sequenceLength: "sequence_length"
      };
      let orderByClause = "ORDER BY created_at DESC";  

      if (sort) {
          const [field, direction] = sort.split(":");
          const dbField = validSortFields[field];
          const sortDirection = (direction && direction.toLowerCase() === "desc") ? "DESC" : "ASC";

          if (dbField) {
              orderByClause = `ORDER BY ${dbField} ${sortDirection}`;
          }
      }

      
      const totalResult = await client.query("SELECT COUNT(*) FROM proteins");
      const total = parseInt(totalResult.rows[0].count, 10);

      
      const result = await client.query(
          `SELECT * FROM proteins ${orderByClause} LIMIT $1 OFFSET $2`,
          [limit, offset]
      );

    
      const proteins = result.rows.map(protein => ({
          proteinId: protein.protein_id,
          name: protein.name,
          description: protein.description,
          molecularWeight: protein.molecular_weight,
          sequenceLength: protein.sequence_length,
          createdAt: protein.created_at,
          updatedAt: protein.updated_at,
          sequenceUrl: protein.sequence_url
      }));

      
      res.status(200).json({
          proteins,
          total,     
          limit,     
          offset     
      });

  } catch (error) {
      console.error("Error fetching proteins:", error);
      res.status(500).json({ error: "Internal Server Error" });
  } finally {
      client.release();
  }
});
//-------------------------------------------------------------------------------------------//

//--------------------------SEARCH----------------------------------------------------------------//
const allowedParams = ["name", "molecularWeight", "sequenceLength", "motif", "sort", "limit", "offset"];
const allowedOperators = ["gt", "gte", "lt", "lte", "eq"];
const allowedSortFields = {
  name: "p.name",
  createdAt: "p.created_at", 
  molecularWeight: "p.molecular_weight",
  sequenceLength: "p.sequence_length",
};

function validateQueryParams(req, res, next) {
    //console.log("Incoming Query Params:", req.query);

    for (const param in req.query) {
        if (!allowedParams.includes(param)) {
            return res.status(400).json({ error: `Invalid query parameter: ${param}` });
        }
    }

    // Validate sorting
    if (req.query.sort) {
        const [field] = req.query.sort.split(":");
        if (!allowedSortFields[field]) {
            return res.status(400).json({ error: `Invalid sort field: ${field}` });
        }
    }

    // Validate molecularWeight
    if (req.query.molecularWeight) {
        for (const op in req.query.molecularWeight) {
            if (!allowedOperators.includes(op) || isNaN(req.query.molecularWeight[op])) {
                return res.status(400).json({ error: `Invalid molecularWeight query` });
            }
        }
    }


    if (req.query.sequenceLength) {
        for (const op in req.query.sequenceLength) {
            if (!allowedOperators.includes(op) || isNaN(req.query.sequenceLength[op])) {
                return res.status(400).json({ error: `Invalid sequenceLength query` });
            }
        }
    }

    next(); 
}

app.get("/api/proteins/search", authenticateUser, validateQueryParams, async (req, res) => {
  const client = await pool.connect();
  try {
      const { name, motif, sort, limit = 10, offset = 0 } = req.query;
      let query = "SELECT DISTINCT p.* FROM proteins p";
      let totalCountQuery = "SELECT COUNT(DISTINCT p.protein_id) FROM proteins p";
      const params = [];
      const conditions = [];

     
      if (req.query.motif) {
        let motif = req.query.motif;
       
        motif = motif.replace(/^\?motif=/, "");
    
       
        const sanitizedMotif = motif.replace(/([.*+?^=!:${}()|[\]/\\])/g, "\\$1");
    
        query += " INNER JOIN fragments f ON p.protein_id = f.protein_id INNER JOIN motifs m ON f.fragment_id = m.fragment_id";
        totalCountQuery += " INNER JOIN fragments f ON p.protein_id = f.protein_id INNER JOIN motifs m ON f.fragment_id = m.fragment_id";
    
        conditions.push(`m.motif_pattern ~ $${params.length + 1}`);

        params.push(motif);

    }
    

      
      if (name) {
          conditions.push(`p.name ILIKE $${params.length + 1}`);
          params.push(`%${name}%`);
      }

      // Apply molecularWeight filtering
      if (req.query.molecularWeight) {
          for (const op in req.query.molecularWeight) {
              const sqlOp = { gt: ">", gte: ">=", lt: "<", lte: "<=", eq: "=" }[op];
              params.push(Number(req.query.molecularWeight[op]));
              conditions.push(`p.molecular_weight ${sqlOp} $${params.length}`);
          }
      }

      if (req.query.sequenceLength) {
          for (const op in req.query.sequenceLength) {
              const sqlOp = { gt: ">", gte: ">=", lt: "<", lte: "<=", eq: "=" }[op];
              params.push(Number(req.query.sequenceLength[op]));
              conditions.push(`p.sequence_length ${sqlOp} $${params.length}`);
          }
      }

      if (conditions.length > 0) {
          query += " WHERE " + conditions.join(" AND ");
          totalCountQuery += " WHERE " + conditions.join(" AND ");
      }

    
      const countParams = [...params];

      if (sort) {
        const [field, direction] = sort.split(":");
    
        
        const dbField = { 
            name: "p.name", 
            createdAt: "p.created_at", 
            molecularWeight: "p.molecular_weight", 
            sequenceLength: "p.sequence_length" 
        }[field];
    
        
        if (!dbField) {
            return res.status(400).json({ error: `Invalid sort field: ${field}` });
        }
    
        const sortDirection = direction?.toLowerCase() === "desc" ? "DESC" : "ASC";
        query += ` ORDER BY ${dbField} ${sortDirection}`;
    }
    

      
      query += ` LIMIT $${params.length + 1} OFFSET $${params.length + 2}`;
      params.push(Number(limit), Number(offset));

      //console.log("Final Query:", query, "Params:", params);
      const { rows } = await client.query(query, params);

      
      const totalCountResult = await client.query(totalCountQuery, countParams);
      const total = parseInt(totalCountResult.rows[0].count, 10);

      
      res.status(200).json({
        message: "200 OK: Successfully retrieved the list of proteins",
        proteins: rows,
        total,
        limit: Number(limit),
        offset: Number(offset),
    });
    

  } catch (error) {
      console.error("Error searching proteins:", error);
      res.status(500).json({ error: "Internal Server Error" });
  } finally {
      client.release();
  }
});
//----------------------------------------------------------------------------------------------------------//



app.get("/api/proteins/:id", async (req, res) => {
  const client = await pool.connect();
  try {
      const { id } = req.params;

      if (!/^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/.test(id)) {
          return res.status(400).json({ error: "Invalid proteinId format" });
      }
    
      const result = await client.query(
          `SELECT * FROM proteins WHERE protein_id = $1`, 
          [id]
      );
    
      if (result.rows.length === 0) {
          return res.status(404).json({ error: "Protein with the given ID does not exist." });
      }
    
      const protein = result.rows[0];   
      res.status(200).json({
          proteinId: protein.protein_id,
          name: protein.name,
          description: protein.description,
          molecularWeight: protein.molecular_weight,
          sequenceLength: protein.sequence_length,
          createdAt: protein.created_at,
          updatedAt: protein.updated_at,
          sequenceUrl: protein.sequence_url
      });

  } catch (error) {
      console.error("Error fetching protein:", error);
      res.status(500).json({ error: "Internal Server Error" });
  } finally {
      client.release();
  }
});


app.get("/api/fragments/:fragmentId", authenticateUser, async (req, res) => {
  const client = await pool.connect();
  try {
    const { fragmentId } = req.params;

   
    if (!/^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/.test(fragmentId)) {
      return res.status(400).json({ error: "Invalid fragmentId format" });
    }

   
    const fragmentQuery = `
      SELECT f.fragment_id, f.protein_id, f.sequence, f.start_position, f.end_position,
             f.secondary_structure, f.created_at, f.url
      FROM fragments f
      WHERE f.fragment_id = $1;
    `;

    const fragmentResult = await client.query(fragmentQuery, [fragmentId]);

  
    if (fragmentResult.rows.length === 0) {
      return res.status(404).json({ error: "Fragment not found" });
    }
    const fragment = fragmentResult.rows[0];

    
    const motifQuery = `
      SELECT mot.motif_pattern, mot.confidence_score
      FROM motifs mot
      WHERE mot.fragment_id = $1
    `;

    const motifResult = await client.query(motifQuery, [fragmentId]);

    
    const motifs = motifResult.rows.map(motif => motif.motif_pattern);
    

    const confidenceScores = motifResult.rows.map(motif => motif.confidence_score);

    const response = {
      fragmentId: fragment.fragment_id,
      proteinId: fragment.protein_id,
      sequence: fragment.sequence,
      startPosition: fragment.start_position,
      endPosition: fragment.end_position,
      motifs,  
      secondaryStructure: fragment.secondary_structure,
      confidenceScores, 
      createdAt: fragment.created_at,
      url: fragment.url
    };

    //res.status(200).json(response);
    res.status(200).json({
      message: "Successfully retrieved the fragment.",
      response
    });

  } catch (error) {
    console.error("Error fetching fragment:", error);
    res.status(500).json({ error: "Internal Server Error" });
  } finally {
    client.release();
  }
});

app.get("/api/proteins/:proteinId/fragments", authenticateUser, async (req, res) => {
  const client = await pool.connect();
  try {
    const { proteinId } = req.params;

    
    if (!/^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/.test(proteinId)) {
      return res.status(400).json({ error: "Invalid proteinId format" });
    }

   
    const proteinCheckQuery = `SELECT protein_id FROM proteins WHERE protein_id = $1`;
    const proteinCheckResult = await client.query(proteinCheckQuery, [proteinId]);

    if (proteinCheckResult.rows.length === 0) {
      return res.status(404).json({ error: "Protein with the given ID does not exist." });
    }

    
    const fragmentsQuery = `
      SELECT f.fragment_id, f.protein_id, f.sequence, f.start_position, f.end_position,
             f.secondary_structure, f.created_at, f.url,
             COALESCE(json_agg(
                json_build_object(
                  'motifPattern', m.motif_pattern, 
                  'confidenceScore', m.confidence_score
                )
             ) FILTER (WHERE m.motif_id IS NOT NULL), '[]') AS motifs
      FROM fragments f
      LEFT JOIN motifs m ON f.fragment_id = m.fragment_id
      WHERE f.protein_id = $1
      GROUP BY f.fragment_id
      ORDER BY f.start_position ASC;
    `;

    const fragmentResults = await client.query(fragmentsQuery, [proteinId]);

   
    const fragments = fragmentResults.rows.map(fragment => ({
      fragmentId: fragment.fragment_id,
      proteinId: fragment.protein_id,
      sequence: fragment.sequence,
      startPosition: fragment.start_position,
      endPosition: fragment.end_position,
      motifs: fragment.motifs.map(motif => motif.motifPattern), 
      secondaryStructure: fragment.secondary_structure,
      confidenceScores: fragment.motifs.map(motif => motif.confidenceScore), 
      createdAt: fragment.created_at,
      url: fragment.url
    }));

    res.status(200).json({
      message: "Successfully retrieved the list of fragments",
      fragments
    });

  } catch (error) {
    console.error("Error fetching fragments for protein:", error);
    res.status(500).json({ error: "Internal Server Error" });
  } finally {
    client.release();
  }
});

app.post("/api/proteins/sequence", authenticateUser, async (req, res) => {
  const client = await pool.connect();

  try {
    const sequence = req.body.trim(); // Extract plain text sequence

   
    if (!sequence || sequence.length > 1000 || !/^[ACDEFGHIKLMNPQRSTVWY]+$/.test(sequence)) {
      return res.status(400).json({
        error: "Invalid input data or sequence length exceeds 1000 characters."
      });
    }

 
    const proteinId = uuidv4();
    const name = `Protein_${sequence.slice(0, 8)}_${Math.floor(Date.now() / 1000)}`;
    const molecularWeight = calculateMolecularWeight(sequence);
    const sequenceLength = sequence.length;
    const createdAt = new Date().toISOString();
    const updatedAt = createdAt;

 
    //const proteinFilePath = `data/proteins/${proteinId}.json`;
    const proteinData = {
      metadata: { version: "1.0", createdAt, updatedAt },
      data: { id: proteinId, name, sequence, description: "", molecularWeight }
    };

    // fs.mkdirSync(path.dirname(proteinFilePath), { recursive: true });
    // fs.writeFileSync(proteinFilePath, JSON.stringify(proteinData, null, 2));

   
    const sequenceUrl = `http://localhost:3000/api/proteins/${proteinId}/sequence`;

    await client.query("BEGIN");

 
    await client.query(
      `INSERT INTO proteins (protein_id, name, description, molecular_weight, sequence_length, sequence_url, created_at, updated_at)
       VALUES ($1, $2, $3, $4, $5, $6, $7, $8)`,
      [proteinId, name, "", molecularWeight, sequenceLength, sequenceUrl, createdAt, updatedAt]
    );

   
    await fragmentAndStoreSequence(client, proteinId, sequence);

    await client.query("COMMIT");

    res.status(201).json({
      message: "Protein created successfully",
      proteinId,
      name,
      molecularWeight,
      sequenceLength,
      createdAt,
      updatedAt,
      sequenceUrl
    });

  } catch (error) {
    await client.query("ROLLBACK");
    console.error("Transaction failed:", error);
    res.status(500).json({ error: "Internal Server Error" });
  } finally {
    client.release();
  }
});

app.get("/api/proteins/:proteinId/sequence", authenticateUser, async (req, res) => {
  const client = await pool.connect();
  try {
    const { proteinId } = req.params;

    
    if (!/^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/.test(proteinId)) {
      return res.status(400).json({ error: "Invalid proteinId format" });
    }

    
    const proteinCheck = await client.query(
      "SELECT name FROM proteins WHERE protein_id = $1",
      [proteinId]
    );
    if (proteinCheck.rows.length === 0) {
      return res.status(404).json({ error: "Protein not found" });
    }

    
    const fragmentResult = await client.query(
      `SELECT sequence FROM fragments 
       WHERE protein_id = $1 
       ORDER BY start_position ASC`,
      [proteinId]
    );

    const fullSequence = fragmentResult.rows.map(row => row.sequence).join("");

    
    const accept = req.get("Accept");

    if (accept === "text/plain") {
      res.type("text/plain").send(fullSequence);
    } else {
      res.json({
        proteinId,
        name: proteinCheck.rows[0].name,
        sequence: fullSequence,
        length: fullSequence.length,
        type: "reconstructed"
      });
    }
  } catch (error) {
    console.error("Error reconstructing protein sequence:", error);
    res.status(500).json({ error: "Internal Server Error" });
  } finally {
    client.release();
  }
});


app.put("/api/proteins/:proteinId", authenticateUser, async (req, res) => {
  const client = await pool.connect();
  try {
      const { proteinId } = req.params;
      const { name, description } = req.body;

      
      if (!/^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/.test(proteinId)) {
          return res.status(400).json({ error: "Invalid proteinId format." });
      }

      
      if (!name && !description) {
          return res.status(400).json({ error: "Invalid input data. Provide at least one updatable field." });
      }
      if (name && name.length > 100) {
          return res.status(400).json({ error: "Protein name exceeds 100 characters." });
      }
      if (description && description.length > 1000) {
          return res.status(400).json({ error: "Description exceeds 1000 characters." });
      }

      
      const checkProtein = await client.query("SELECT * FROM proteins WHERE protein_id = $1", [proteinId]);
      if (checkProtein.rows.length === 0) {
          return res.status(404).json({ error: "Protein with given ID does not exist." });
      }

      
      const updates = [];
      const params = [];
      if (name) {
          updates.push("name = $1");
          params.push(name);
      }
      if (description) {
          updates.push("description = $" + (params.length + 1));
          params.push(description);
      }

      params.push(new Date().toISOString(), proteinId);     
      const updateQuery = `
          UPDATE proteins
          SET ${updates.join(", ")}, updated_at = $${params.length - 1}
          WHERE protein_id = $${params.length}
          RETURNING *;
      `;
      const result = await client.query(updateQuery, params);

      res.status(200).json({
          message: "Protein updated successfully",
          protein: result.rows[0],
      });

  } catch (error) {
      console.error("Error updating protein:", error);
      res.status(500).json({ error: "Internal Server Error" });
  } finally {
      client.release();
  }
});

app.delete("/api/proteins/:proteinId", authenticateUser, async (req, res) => {
  const client = await pool.connect();
  try {
      const { proteinId } = req.params;

    
      if (!/^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/.test(proteinId)) {
          return res.status(400).json({ error: "Invalid proteinId format" });
      }

      await client.query("BEGIN");

     
      const proteinCheck = await client.query("SELECT * FROM proteins WHERE protein_id = $1", [proteinId]);

      if (proteinCheck.rows.length === 0) {
          await client.query("ROLLBACK");
          return res.status(404).json({ error: "Protein with given ID does not exist" });
      }

      
      await client.query("DELETE FROM proteins WHERE protein_id = $1", [proteinId]);

      await client.query("COMMIT");

      res.status(204).send(); 
  } catch (error) {
      await client.query("ROLLBACK");
      console.error("Error deleting protein:", error);
      res.status(500).json({ error: "Internal Server Error" });
  } finally {
      client.release();
  }
});


app.get("/api/proteins/:proteinId/structure", authenticateUser, async (req, res, next) => {
  const client = await pool.connect();

  try {
      const { proteinId } = req.params;

   
      if (!/^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/.test(proteinId)) {
          return res.status(400).json({ error: "Invalid proteinId format" });
      }

  
      const proteinCheck = await client.query(
          "SELECT name FROM proteins WHERE protein_id = $1",
          [proteinId]
      );
      if (proteinCheck.rows.length === 0) {
          return res.status(404).json({ error: "Protein not found" });
      }

    
      const fragmentResult = await client.query(
          `SELECT sequence FROM fragments 
           WHERE protein_id = $1 
           ORDER BY start_position ASC`,
          [proteinId]
      );

      const sequence = fragmentResult.rows.map(row => row.sequence).join("");

      if (!sequence || sequence.length === 0) {
          return res.status(404).json({ error: "No sequence data found for this protein" });
      }

   
      const structure = predictSecondaryStructure(proteinId, sequence);

      if (req.accepts("application/json")) {
          return res.json(structure);
      } else if (req.accepts("image/svg+xml")) {
          const svg = generateStructureSVG(sequence, structure.secondaryStructure);
          return res.type("image/svg+xml").send(svg);
      } else {
          return res.status(406).json({ error: "406 Not Acceptable: Requested content type not available" });
      }

  } catch (error) {
      console.error("Error in /structure endpoint:", error);
      next(error);
  } finally {
      client.release();
  }
});




function generateStructureSVG(sequence, secondaryStructure) {
  const svgWidth = sequence.length * 10;  
  const svgHeight = 50;  

  let svg = `<svg width="${svgWidth}" height="${svgHeight}" xmlns="http://www.w3.org/2000/svg">`;

 
  for (let i = 0; i < sequence.length; i++) {
      let color;
      switch (secondaryStructure[i]) {
          case 'H': color = 'red'; break;    
          case 'E': color = 'yellow'; break;  
          case 'C': color = 'gray'; break;    
          default: color = 'black';            
      }

      svg += `<rect x="${i * 10}" y="0" width="10" height="30" fill="${color}"/>`;
  }

  // Add legend
      svg += `
  <rect x="10" y="35" width="10" height="10" fill="red" />
  <text x="25" y="45" font-size="10">alpha-helix</text>
  <rect x="70" y="35" width="10" height="10" fill="yellow" />
  <text x="85" y="45" font-size="10">beta-strand</text>
  <rect x="140" y="35" width="10" height="10" fill="gray" />
  <text x="155" y="45" font-size="10">coil</text>
  `;

  svg += `</svg>`;
  return svg;
}


app.use((err, req, res, next) => {
  if (err instanceof SyntaxError && err.status === 400 && "body" in err) {
      return res.status(400).json({
          error: "Invalid JSON format",
          details: { message: "Ensure the request body is a properly formatted JSON object." }
      });
  }
  next();
});


function errorHandler(err, req, res, next) {
  console.error(err);
  if (err instanceof NotFoundError) {
      res.status(404).json({ error: err.message });
  } else if (err instanceof ConflictError) {
      res.status(409).json({ error: err.message });
  } else {
      res.status(500).json({ error: "Internal server error" });
  }
}

app.use(errorHandler);




