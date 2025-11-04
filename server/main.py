from fastapi import FastAPI, File, UploadFile, Form, Request
from fastapi.responses import Response, JSONResponse
from fastapi.templating import Jinja2Templates
from contextlib import asynccontextmanager
import csv
from starlette.responses import FileResponse
from io import StringIO
from Bio import SeqIO
from choppy.homology_finder import process_background_sequences, process_query_sequences, find_non_homologous_regions, create_annotated_record, load_trie
from choppy.fragment_annotator import annotate_fragments, FragmentConfig
import base64
from pathlib import Path

trie_cache = {}
kmer_sizes = {}

@asynccontextmanager
async def lifespan(app: FastAPI):
    with open("data/data_list.csv", "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            full_path = "data/" + row["file_name"]
            if Path(full_path).exists():
                kmer_sizes[row["id"]] = int(row["kmer_size"])
                trie_cache[row["id"]] = load_trie(full_path)

    yield
    trie_cache.clear()
    kmer_sizes.clear()

app = FastAPI(root_path="/choppy", lifespan=lifespan)
templates = Jinja2Templates(directory="server/static")
@app.get("/")
async def read_root(request: Request):
    with open("data/data_list.csv", "r") as f:
        reader = csv.DictReader(f)
        bgs = [row for row in reader]
    return templates.TemplateResponse("index.html", {"request": request, "bgs": bgs})

@app.post("/process-and-fragment")
async def process_and_fragment_endpoint(
    query_file: UploadFile = File(...),
    custom_background_files: list[UploadFile] = File(default=[]),
    min_size: int = Form(500),
    max_size: int = Form(3000),
    min_overlap: int = Form(60),
    max_overlap: int = Form(100),
    boundary_motif: str = Form(""),
    min_step: int = Form(10),
    precalculated_bgs: str = Form("")
):
    """
    Combined endpoint that:
    1. Annotates homology-free regions
    2. Fragments the sequence based on those regions
    Returns both the annotated GenBank file and fragments FASTA file as JSON
    """
    print(f"Query file: {query_file.filename}")
    print(f"Background files: {[f.filename for f in custom_background_files]}")
    print(f"Fragment params - Size: {min_size}-{max_size}, Overlap: {min_overlap}-{max_overlap}")
    print(f"Boundary motif: {boundary_motif or 'None'}, Min step: {min_step}")
    print(f"Precalculated backgrounds: {precalculated_bgs}")

    #temporary
    kmer_size = 15
    threshold = 60
    try:
        query_content = await query_file.read()
        query_text = query_content.decode('utf-8')
        query_io = StringIO(query_text)
        
        query_format = "genbank" if query_file.filename.endswith((".gb", ".gbk", ".genbank")) else "fasta"
        query_sequences = list(SeqIO.parse(query_io, query_format))
        
        background_sequences = []
        for bg_file in custom_background_files:
            bg_content = await bg_file.read()
            format_type = "fasta" if bg_file.filename.endswith((".fa", ".fasta")) else "genbank"
            bg_text = bg_content.decode('utf-8')
            bg_io = StringIO(bg_text)
            background_sequences.extend(list(SeqIO.parse(bg_io, format_type)))
                
        bg_tries = process_background_sequences(background_sequences, kmer_size, merge=False)
        for bg_id in precalculated_bgs.split(","):
            if bg_id in trie_cache:
                bg_tries.append(trie_cache[bg_id])
        
        query_tries = process_query_sequences(query_sequences, kmer_size)

        print(f"Processing {len(query_sequences)} query sequence(s) with {len(bg_tries)} background sequence(s)...")
        annotated_records = []
        
        for query_seq in query_sequences:
            regions = find_non_homologous_regions(
                query_seq, query_tries[query_seq.id], bg_tries, kmer_size, threshold
            )
            annotated_record = create_annotated_record(query_seq, regions)
            annotated_records.append(annotated_record)
        
        print("Homology annotation complete!")
        
        # For the time being, it is assumed that only one query sequence is processed
        print("Running fragmentor (placeholder)...")
        
        config = FragmentConfig(
            min_size=min_size,
            max_size=max_size,
            min_overlap=min_overlap,
            max_overlap=max_overlap,
            motif=boundary_motif if boundary_motif else None,
            min_step=min_step
        )
        records, frag_records = annotate_fragments(annotated_records[0], config)
        
        # Write fragments to FASTA StringIO
        fragments_io = StringIO()
        SeqIO.write(frag_records, fragments_io, "fasta")
        fragments_content = fragments_io.getvalue()

        # Write annotated GenBank file
        annotated_io = StringIO()
        SeqIO.write(records, annotated_io, "genbank")
        annotated_content = annotated_io.getvalue()


        print("Fragmentation complete (placeholder)!")
        
        annotated_b64 = base64.b64encode(annotated_content.encode('utf-8')).decode('utf-8')
        fragments_b64 = base64.b64encode(fragments_content.encode('utf-8')).decode('utf-8')
        
        return JSONResponse({
            "status": "success",
            "annotated_file": annotated_b64,
            "annotated_filename": "annotated_sequences.gb",
            "fragments_file": fragments_b64,
            "fragments_filename": "fragments.fasta",
            "message": "Processing complete!"
        })
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return JSONResponse({
            "status": "error",
            "message": str(e)
        }, status_code=500)
