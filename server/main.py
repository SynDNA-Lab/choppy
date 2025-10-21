from fastapi import FastAPI, File, UploadFile, Form
from fastapi.responses import Response
from starlette.responses import FileResponse
from io import BytesIO, StringIO
from Bio import SeqIO
from tarnche.homology_finder import process_background_sequences, process_query_sequences, find_non_homologous_regions, create_annotated_record


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