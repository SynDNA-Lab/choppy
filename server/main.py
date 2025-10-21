from fastapi import FastAPI, File, UploadFile, Form
from fastapi.responses import Response
from starlette.responses import FileResponse
from io import BytesIO
from Bio import SeqIO
from tarnche.homology_finder import annotate_homology


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

    query_content = await query_file.read()
    query_io = BytesIO(query_content)
    
    background_ios = []
    for bg_file in background_files:
        bg_content = await bg_file.read()
        background_ios.append(BytesIO(bg_content))

    print(f"Processing with {len(background_ios)} background sequences...")

    annotated_records = annotate_homology(
        query_io,
        background_ios,
        kmer_size=kmer_size,
        threshold=threshold,
        output_path="temp_output.gb"
    )

    # print(f"Annotation complete, preparing response...")

    # output_io = BytesIO()
    # SeqIO.write(annotated_records, output_io, "genbank")
    # output_content = output_io.getvalue()

    # return Response(
    #     content=output_content,
    #     media_type="application/octet-stream",
    #     headers={"Content-Disposition": "attachment; filename=annotated_sequences.gb"}
    # )