// main.nf
// Nextflow Pipeline for Phage Isolate Assembly & Annotation
// Workflow: QC -> Flye -> Medaka -> CheckV/ViralVerify -> Pharokka -> Phold -> Phynteny -> Report
// Usage:
//   Real run:  nextflow run main.nf -c nextflow.config --inputFile long_reads.fastq 
//   Test run:  nextflow run main.nf -c nextflow.config --inputFile long_reads.fastq -stub

// --- Parameters ---
params.inputFile  = null
params.outdir     = './results'

// ============================================================
// 1. QC & Trimming: fastplong
// ============================================================
process QC_TRIM {
    publishDir "${params.outdir}/01_qc", mode: 'copy'

    input:
    path reads

    output:
    path "*_trimmed.fastq.gz", emit: trimmed_reads
    path "fastplong.html",  emit: qc_report
    path "fastplong.json",  emit: qc_report_json

    script:
    """
    fastplong \
        -i ${reads} \
        -o ${reads.baseName}_trimmed.fastq.gz \
        --n_base_limit 500 \
        --thread ${task.cpus}
    """

    stub:
    """
    echo "Dummy trimmed content" | gzip > ${reads.baseName}_trimmed.fastq.gz
    echo "<html><body>Stub QC Report</body></html>" > fastplong.html
    echo '{"summary": "stub"}' > fastplong.json
    """
}

// ============================================================
// 2. Assembly: Flye 
// ============================================================
process ASSEMBLY {
    publishDir "${params.outdir}/02_assembly", mode: 'copy'

    input:
    path trimmed_reads

    output:
    path "flye_out/assembly.fasta",       emit: assembly
    path "flye_out/assembly_info.txt",    emit: assembly_info

    script:
    """
    flye \
        --nano-corr ${trimmed_reads} \
        --out-dir flye_out \
        --meta \
        --threads ${task.cpus}
    """

    stub:
    """
    mkdir -p flye_out
    printf '>contig_1 length=50000\\nATCGATCGATCG\\n' > flye_out/assembly.fasta
    echo "assembly_info: stubbed" > flye_out/assembly_info.txt
    """
}

// ============================================================
// 3. Polishing: Medaka
// ============================================================
process POLISH {
    publishDir "${params.outdir}/02_assembly", mode: 'copy'

    input:
    path assembly
    path trimmed_reads

    output:
    path "medaka_out/consensus.fasta", emit: polished_assembly

    script:
    """
    medaka_consensus \
        -i ${trimmed_reads} \
        -d ${assembly} \
        -o medaka_out \
        -t ${task.cpus}
    """

    stub:
    """
    mkdir -p medaka_out
    printf '>contig_1_polished length=50000\\nATCGATCGATCG\\n' > medaka_out/consensus.fasta
    """
}

// ============================================================
// 4. Assembly QC: CheckV
// ============================================================
process CHECKV {
    publishDir "${params.outdir}/03_qc", mode: 'copy'

    input:
    path polished_assembly

    output:
    path "checkv_out", emit: checkv_results

    script:
    """
    checkv end_to_end ${polished_assembly} checkv_out -t ${task.cpus} -d ~/micromamba/envs/phage-pipeline/db/checkv-db-v1.5
    """

    stub:
    """
    mkdir -p checkv_out
    printf 'contig_id,completeness,contamination\\nviral_contig_1,95.5,1.2\\n' > checkv_out/quality_summary.tsv
    """
}

// ============================================================
// 5. Taxonomic Verification: ViralVerify
// ============================================================
process VIRAL_VERIFY {
    publishDir "${params.outdir}/03_qc", mode: 'copy'

    input:
    path polished_assembly

    output:
    path "viralverify_out", emit: viralverify_dir

    script:
    """
    viralverify -f ${polished_assembly} -o viralverify_out --hmm ~/micromamba/envs/phage-pipeline/db/nbc_hmms.hmm
    """

    stub:
    """
    printf 'Contig\\tPredicted_Label\\tFamily\\tGenus\\tviral_score\\n' > viralverify_report.txt
    printf 'contig_1\\tVirus\\tSiphoviridae\\tUnknown\\t0.98\\n' >> viralverify_report.txt
    """
}

// ============================================================
// 6a. Phage Annotation: Pharokka
// ============================================================
process PHAROKKA {
    publishDir "${params.outdir}/04_annotation", mode: 'copy'

    input:
    path polished_assembly

    output:
    path "pharokka_out/pharokka.gff",         emit: pharokka_gff
    path "pharokka_out/pharokka.gbk",         emit: pharokka_gbk
    path "pharokka_out",                      emit: pharokka_dir

    script:
    """
    pharokka.py \
        -i ${polished_assembly} \
        -o pharokka_out \
        -t ${task.cpus} \
        -d ~/micromamba/envs/phage-pipeline/db/ \
        -f
    """

    stub:
    """
    mkdir -p pharokka_out
    printf '##gff-version 3\\nviral_contig_1\\tpharokka\\tCDS\\t1\\t500\\t.\\t+\\t0\\tID=CDS_1;phrog=phrog_1;function=tail fiber\\n' \
        > pharokka_out/pharokka.gff
    echo "LOCUS viral_contig_1" > pharokka_out/pharokka.gbk
    """
}

// // ============================================================
// // 6b. Structure-based Annotation: Phold
// // ============================================================
// process PHOLD {
//     publishDir "${params.outdir}/04_annotation/phold", mode: 'copy'

//     input:
//     path gbk_file

//     output:
//     path "phold_out/phold.gff",         emit: phold_gff
//     path "phold_out/phold.gbk",         emit: phold_gbk
//     path "phold_out",                   emit: phold_dir

//     script:
//     """
//     phold run \
//         -i ${gbk_file} \
//         -o phold_out \
//         -t ${task.cpus} \
//         -f
//     """

//     stub:
//     """
//     mkdir -p phold_out
//     printf '##gff-version 3\\nviral_contig_1\\tphold\\tCDS\\t1\\t500\\t.\\t+\\t0\\tID=CDS_1;function=tail fiber protein (3D)\\n' \
//         > phold_out/phold.gff
//     cp ${gbk_file} phold_out/phold.gbk
//     """
// }

// // ============================================================
// // 6c. Synteny-based Context: Phynteny
// // ============================================================
// process PHYNTENY {
//     publishDir "${params.outdir}/04_annotation/phynteny", mode: 'copy'

//     input:
//     path gbk_file

//     output:
//     path "phynteny_out/phynteny.gbk",         emit: phynteny_gbk
//     path "phynteny_out/phynteny_summary.tsv", emit: phynteny_summary

//     script:
//     """
//     mkdir -p phynteny_input
//     cp ${gbk_file} phynteny_input/
    
//     phynteny \
//         phynteny_input \
//         phynteny_out
//     """

//     stub:
//     """
//     mkdir -p phynteny_out
//     echo "LOCUS viral_contig_1_phynteny" > phynteny_out/phynteny.gbk
//     printf 'contig\\tsynteny_group\\trelated_phages\\ncontig_1\\tSG_042\\t12\\n' \
//         > phynteny_out/phynteny_summary.tsv
//     """
// }

// // ============================================================
// // 7. Generate PDF Report
// // ============================================================
// process GENERATE_REPORT {
//     publishDir "${params.outdir}/05_report", mode: 'copy'

//     input:
//     path qc_json
//     path assembly_info
//     path checkv_results
//     path viralverify_dir
//     path pharokka_dir

//     output:
//     path "phage_analysis_report.pdf", emit: pdf_report

//     script:
//     """
//     #!/usr/bin/env python3
//     from reportlab.lib.pagesizes import letter
//     from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
//     from reportlab.lib.units import inch
//     from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
//     from reportlab.lib import colors
//     from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
//     import json
//     import os
//     from datetime import datetime

//     # Read input data
//     qc_data = json.load(open("${qc_json}"))
    
//     # Parse CheckV results
//     checkv_summary = None
//     for f in os.listdir("${checkv_results}"):
//         if "quality_summary" in f:
//             checkv_summary = f
//             break
    
//     # Parse Pharokka summary
//     pharokka_summary = None
//     for f in os.listdir("${pharokka_dir}"):
//         if f.endswith("_cds_functions.tsv"):
//             pharokka_summary = f
//             break

//     # Create PDF
//     doc = SimpleDocTemplate("phage_analysis_report.pdf", pagesize=letter,
//                            rightMargin=72, leftMargin=72,
//                            topMargin=72, bottomMargin=18)
    
//     story = []
//     styles = getSampleStyleSheet()
    
//     # Custom styles
//     title_style = ParagraphStyle(
//         'CustomTitle',
//         parent=styles['Heading1'],
//         fontSize=24,
//         textColor=colors.HexColor('#1a1a1a'),
//         spaceAfter=30,
//         alignment=TA_CENTER
//     )
    
//     heading_style = ParagraphStyle(
//         'CustomHeading',
//         parent=styles['Heading2'],
//         fontSize=14,
//         textColor=colors.HexColor('#2c5aa0'),
//         spaceAfter=12,
//         spaceBefore=12
//     )
    
//     # Title
//     story.append(Paragraph("Phage Isolate Analysis Report", title_style))
//     story.append(Spacer(1, 0.2*inch))
    
//     # Metadata
//     date_str = datetime.now().strftime("%Y-%m-%d %H:%M")
//     story.append(Paragraph(f"<b>Date:</b> {date_str}", styles['Normal']))
//     story.append(Paragraph(f"<b>Input File:</b> ${params.inputFile}", styles['Normal']))
//     story.append(Spacer(1, 0.3*inch))
    
//     # Methods Section
//     story.append(Paragraph("Methods", heading_style))
//     methods_text = \"\"\"
//     <b>Quality Control:</b> Raw Oxford Nanopore reads were processed using fastplong 
//     to remove adapters and low-quality bases (N-base limit: 500).
//     <br/><br/>
//     <b>Assembly:</b> Trimmed reads were assembled using Flye v2.9 in metagenomic mode 
//     (--meta) with the --nano-corr preset for error-corrected reads.
//     <br/><br/>
//     <b>Polishing:</b> The draft assembly was polished using Medaka consensus calling 
//     to correct remaining errors in homopolymer regions and improve base accuracy.
//     <br/><br/>
//     <b>Quality Assessment:</b> Viral completeness and contamination were assessed using 
//     CheckV v1.5. Taxonomic verification was performed with viralVerify to confirm 
//     viral origin and predict taxonomic classification.
//     <br/><br/>
//     <b>Functional Annotation:</b> Gene calling and functional annotation were performed 
//     using Pharokka, which identifies coding sequences and assigns functions based on 
//     the PHROGs database of phage orthologous groups.
//     \"\"\"
//     story.append(Paragraph(methods_text, styles['Normal']))
//     story.append(PageBreak())
    
//     # Results Section
//     story.append(Paragraph("Results", heading_style))
    
//     # QC Stats
//     story.append(Paragraph("<b>1. Quality Control Summary</b>", styles['Heading3']))
//     qc_table_data = [
//         ["Metric", "Value"],
//         ["Total Reads", str(qc_data.get("total_reads", "N/A"))],
//         ["Total Bases", str(qc_data.get("total_bases", "N/A"))],
//         ["Mean Read Length", str(qc_data.get("mean_length", "N/A"))],
//         ["Mean Quality", str(qc_data.get("mean_quality", "N/A"))]
//     ]
//     qc_table = Table(qc_table_data, colWidths=[3*inch, 2*inch])
//     qc_table.setStyle(TableStyle([
//         ('BACKGROUND', (0,0), (-1,0), colors.grey),
//         ('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke),
//         ('ALIGN', (0,0), (-1,-1), 'LEFT'),
//         ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
//         ('FONTSIZE', (0,0), (-1,0), 12),
//         ('BOTTOMPADDING', (0,0), (-1,0), 12),
//         ('BACKGROUND', (0,1), (-1,-1), colors.beige),
//         ('GRID', (0,0), (-1,-1), 1, colors.black)
//     ]))
//     story.append(qc_table)
//     story.append(Spacer(1, 0.3*inch))
    
//     # Assembly Stats
//     story.append(Paragraph("<b>2. Assembly Statistics</b>", styles['Heading3']))
//     assembly_text = "Assembly completed successfully. See assembly_info.txt for contig-level details."
//     story.append(Paragraph(assembly_text, styles['Normal']))
//     story.append(Spacer(1, 0.2*inch))
    
//     # CheckV Results
//     if checkv_summary:
//         story.append(Paragraph("<b>3. Viral Quality Assessment (CheckV)</b>", styles['Heading3']))
//         checkv_path = os.path.join("${checkv_results}", checkv_summary)
//         with open(checkv_path) as f:
//             lines = f.readlines()
//             if len(lines) > 1:
//                 header = lines[0].strip().split('\\t')[:5]  # First 5 columns
//                 data = lines[1].strip().split('\\t')[:5]
//                 checkv_table_data = [header, data]
//                 checkv_table = Table(checkv_table_data)
//                 checkv_table.setStyle(TableStyle([
//                     ('BACKGROUND', (0,0), (-1,0), colors.grey),
//                     ('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke),
//                     ('ALIGN', (0,0), (-1,-1), 'LEFT'),
//                     ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
//                     ('FONTSIZE', (0,0), (-1,-1), 10),
//                     ('GRID', (0,0), (-1,-1), 1, colors.black)
//                 ]))
//                 story.append(checkv_table)
//     story.append(Spacer(1, 0.3*inch))
    
//     # Pharokka Results
//     if pharokka_summary:
//         story.append(Paragraph("<b>4. Functional Annotation (Pharokka)</b>", styles['Heading3']))
//         pharokka_path = os.path.join("${pharokka_dir}", pharokka_summary)
//         with open(pharokka_path) as f:
//             total_cds = sum(1 for line in f) - 1  # minus header
//         story.append(Paragraph(f"Total CDS predicted: {total_cds}", styles['Normal']))
//         story.append(Paragraph("See pharokka output directory for detailed gene annotations and functional categories.", styles['Normal']))
    
//     # Build PDF
//     doc.build(story)
//     """

//     stub:
//     """
//     #!/usr/bin/env python3
//     from reportlab.lib.pagesizes import letter
//     from reportlab.platypus import SimpleDocTemplate, Paragraph
//     from reportlab.lib.styles import getSampleStyleSheet
    
//     doc = SimpleDocTemplate("phage_analysis_report.pdf", pagesize=letter)
//     story = []
//     styles = getSampleStyleSheet()
//     story.append(Paragraph("STUB: Phage Analysis Report", styles['Title']))
//     story.append(Paragraph("This is a stub report for testing.", styles['Normal']))
//     doc.build(story)
//     """
// }

// ============================================================
// Workflow
// ============================================================
workflow {

    if (!params.inputFile) {
        error "Please provide an input file: --inputFile <path_to_fastq>"
    }

    reads_ch = Channel.fromPath(params.inputFile, checkIfExists: true)

    // 1. QC
    qc_out       = QC_TRIM(reads_ch)

    // 2. Assembly
    asm_out      = ASSEMBLY(qc_out.trimmed_reads)

    // 3. Polishing
    polish_out   = POLISH(asm_out.assembly, qc_out.trimmed_reads)

    // 4. QC & Verification (Parallel)
    checkv_out   = CHECKV(polish_out.polished_assembly)
    verify_out   = VIRAL_VERIFY(polish_out.polished_assembly)

    // 5. Annotation Chain
    pharokka_out = PHAROKKA(polish_out.polished_assembly)
    // phold_out    = PHOLD(pharokka_out.pharokka_gbk)
    // phynteny_out = PHYNTENY(phold_out.phold_gbk)

    // // 6. Generate Report
    // report_out = GENERATE_REPORT(
    //     qc_out.qc_report_json,
    //     asm_out.assembly_info,
    //     checkv_out.checkv_results,
    //     verify_out.viralverify_dir,
    //     pharokka_out.pharokka_dir
    // )

    // report_out.pdf_report.view { path ->
    //     """
    //     =========================================
    //     Pipeline completed!
    //     PDF Report: ${path}
    //     =========================================
    //     """
    // }
    
}