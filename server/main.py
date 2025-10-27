from fastapi import FastAPI, File, UploadFile, Form
from fastapi.responses import Response, JSONResponse
from starlette.responses import FileResponse
from io import BytesIO, StringIO
from Bio import SeqIO
from choppy.homology_finder import process_background_sequences, process_query_sequences, find_non_homologous_regions, create_annotated_record
from choppy.fragment_annotator import annotate_fragments, FragmentConfig
import base64


app = FastAPI()

@app.get("/")
async def read_root():
    return FileResponse("server/static/index.html")

@app.post("/annotate-homology")
async def annotate_homology_endpoint(
    query_file: UploadFile = File(...),
    background_files: list[UploadFile] = File(default=[]),
    kmer_size: int = Form(20),
    threshold: int = Form(60),
):
    print(f"Query file: {query_file.filename}")
    print(f"Background files: {[f.filename for f in background_files]}")
    print(f"K-mer size: {kmer_size}, Threshold: {threshold}")

    # Read query file - decode to text for GenBank format
    query_content = await query_file.read()
    query_text = query_content.decode('utf-8')
    query_io = StringIO(query_text)
    query_sequences = list(SeqIO.parse(query_io, "genbank"))

    background_ios = []
    background_sequences = []
    for bg_file in background_files:
        bg_content = await bg_file.read()
        format_type = "fasta" if bg_file.filename.endswith((".fa", ".fasta")) else "genbank"
        bg_text = bg_content.decode('utf-8')
        bg_io = StringIO(bg_text)        
        background_sequences.extend(list(SeqIO.parse(bg_io, format_type)))

    print(f"Processing with {len(background_ios)} background sequences...")

    bg_trie = process_background_sequences(background_sequences, kmer_size)
    query_tries = process_query_sequences(query_sequences, kmer_size)
    output_records = []
    for query_seq in query_sequences:
        regions = find_non_homologous_regions(
            query_seq, query_tries[query_seq.id], bg_trie, kmer_size, threshold
        )
        output_record = create_annotated_record(query_seq, regions)
        output_records.append(output_record)

    output_format = "genbank"
    filename = "annotated_sequences.gb"

    print(f"Annotation complete, preparing response...")

    output_io = StringIO()
    SeqIO.write(output_records, output_io, output_format)
    output_content = output_io.getvalue()

    return Response(
        content=output_content,
        media_type="application/octet-stream",
        headers={"Content-Disposition": f"attachment; filename={filename}"}
    )

@app.post("/process-and-fragment")
async def process_and_fragment_endpoint(
    query_file: UploadFile = File(...),
    background_files: list[UploadFile] = File(default=[]),
    kmer_size: int = Form(20),
    threshold: int = Form(60),
    min_size: int = Form(500),
    max_size: int = Form(3000),
    min_overlap: int = Form(60),
    max_overlap: int = Form(100),
    boundary_motif: str = Form(""),
    min_step: int = Form(10),
):
    """
    Combined endpoint that:
    1. Annotates homology-free regions
    2. Fragments the sequence based on those regions
    Returns both the annotated GenBank file and fragments FASTA file as JSON
    """
    print(f"Query file: {query_file.filename}")
    print(f"Background files: {[f.filename for f in background_files]}")
    print(f"Homology params - K-mer: {kmer_size}, Threshold: {threshold}")
    print(f"Fragment params - Size: {min_size}-{max_size}, Overlap: {min_overlap}-{max_overlap}")
    print(f"Boundary motif: {boundary_motif or 'None'}, Min step: {min_step}")
    
    try:
        query_content = await query_file.read()
        query_text = query_content.decode('utf-8')
        query_io = StringIO(query_text)
        
        query_format = "genbank" if query_file.filename.endswith((".gb", ".gbk", ".genbank")) else "fasta"
        query_sequences = list(SeqIO.parse(query_io, query_format))
        
        background_sequences = []
        for bg_file in background_files:
            bg_content = await bg_file.read()
            format_type = "fasta" if bg_file.filename.endswith((".fa", ".fasta")) else "genbank"
            bg_text = bg_content.decode('utf-8')
            bg_io = StringIO(bg_text)
            background_sequences.extend(list(SeqIO.parse(bg_io, format_type)))
        
        print(f"Processing {len(query_sequences)} query sequence(s) with {len(background_sequences)} background sequence(s)...")
        
        bg_trie = process_background_sequences(background_sequences, kmer_size)
        query_tries = process_query_sequences(query_sequences, kmer_size)
        annotated_records = []
        
        for query_seq in query_sequences:
            regions = find_non_homologous_regions(
                query_seq, query_tries[query_seq.id], bg_trie, kmer_size, threshold
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
