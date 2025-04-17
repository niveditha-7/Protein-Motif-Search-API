# Protein-Motif-Search-API
# Protein Analysis & Motif Detection API

A Dockerized Node.js-based REST API for analyzing and storing protein sequences. It performs metadata extraction (molecular weight, sequence length), fragments sequences, detects motifs, predicts secondary structure, and stores results in a PostgreSQL database.

---

## Features

- Submit protein sequences and store metadata
- Predict secondary structure using GOR method
- Fragment sequences with overlap
- Detect biological motifs:
  - N-glycosylation site: `N[^P][ST][^P]`
  - Casein kinase II phosphorylation site: `[ST].{2}[DE]`
  - Tyrosine kinase phosphorylation site: `[RK].{0,2}[DE]`
- Store and retrieve protein and fragment data from PostgreSQL
- Sequence stored in filesystem, not database
- Fully containerized using Docker

---

## Technologies Used

- Node.js (Express.js)
- PostgreSQL (via Docker)
- Docker for containerization
- REST API + JSON
- Regex for motif matching

---

## Project Structure

```
my_project/
├── node.js               # Main server file
├── .env                  # Environment configuration
├── public/
│   └── sequences/        # Stores sequence text files
├── data/                 # Legacy JSON storage (not used anymore)
├── Dockerfile            # (Optional) If Dockerizing Node app
└── README.md
```

---

## Setup Instructions

### 1. Clone the repository
```bash
git clone <repo-url>
cd "enter into project folder"
```

### 2. Create `.env` file
```env
PORT=3000
MAX_PROTEIN_LENGTH=2000
PG_HOST=localhost
PG_PORT=5432
PG_DATABASE=protein_db
PG_USER=postgres
PG_PASSWORD=password
```

### 3. Run PostgreSQL in Docker
```bash
docker run -d \
  --name postgres \
  --network host \
  -e POSTGRES_PASSWORD=password \
  -e POSTGRES_USER=postgres \
  -e POSTGRES_DB=protein_db \
  -v pgdata:/var/lib/postgresql/data \
  postgres:latest
```

### 4. Connect to the DB and run schema
```bash
docker exec -it postgres psql -U postgres -d protein_db
```
Run the following SQL (or load from schema.sql):
```sql
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";

CREATE TABLE proteins (
  protein_id UUID DEFAULT uuid_generate_v4() PRIMARY KEY,
  name VARCHAR(100) NOT NULL,
  description VARCHAR(1000),
  molecular_weight FLOAT CHECK (molecular_weight > 0),
  sequence_length INTEGER,
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  sequence_url VARCHAR(255)
);

CREATE TABLE fragments (
  fragment_id UUID DEFAULT uuid_generate_v4() PRIMARY KEY,
  protein_id UUID REFERENCES proteins(protein_id) ON DELETE CASCADE,
  sequence VARCHAR(50) CHECK (sequence ~ '^[A-Z]{2,50}$'),
  start_position INTEGER,
  end_position INTEGER,
  secondary_structure VARCHAR(50) CHECK (secondary_structure ~ '^[HEC]+$'),
  confidence_scores JSONB,
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  url VARCHAR(255)
);

CREATE TABLE motifs (
  motif_id UUID DEFAULT uuid_generate_v4() PRIMARY KEY,
  fragment_id UUID REFERENCES fragments(fragment_id) ON DELETE CASCADE,
  motif_pattern VARCHAR(50) NOT NULL,
  motif_type VARCHAR(50),
  start_position INTEGER,
  end_position INTEGER,
  confidence_score FLOAT CHECK (confidence_score >= 0 AND confidence_score <= 1),
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

---

## Running the Node Server
```bash
npm install
node node.js
```

---

## API Endpoints

### POST `/api/proteins`
Create new protein
```json
{
  "sequence": "ACDEFGHIKLMNPQRSTVWY",
  "name": "Test Protein",
  "description": "Example"
}
```
Returns: metadata + sequence URL + fragments stored

---

### GET `/api/proteins`
Paginated list of proteins
- Query params: `limit`, `offset`

---

### GET `/api/proteins/:proteinId`
Return one protein by ID

---

### GET `/api/proteins/:proteinId/fragments`
Get all fragments for a protein

---

### GET `/api/fragments/:fragmentId`
Return fragment + structure + motifs

---

### DELETE `/api/proteins/:proteinId`
Delete protein and all associated fragments and motifs

---

## Testing with Postman or cURL
```bash
curl -X POST http://localhost:3000/api/proteins \
     -H "Content-Type: application/json" \
     -d '{
           "sequence": "ACDEFGHIKLMNPQRSTVWY",
           "name": "My Protein",
           "description": "For testing"
         }'
```

---

## Author & Credits
Developed as part of a bioinformatics backend system coursework.
Special thanks to PostgreSQL documentation, Docker community, and Node.js community.

