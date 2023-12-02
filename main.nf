#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// All of the default parameters are being set in `nextflow.config`

// Import sub-workflows
include { build_solvents } from './modules/system_builder'
include { train_model } from './modules/model_builder'


// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run [-profile local/docker/singularity/slurm] . [--database_path /Path/To/CSV --num_points num_points_coex_curve]

# FreeSolv
nextflow run . -profile docker

# User defined database
nextflow run . -profile docker --database_path databases/spce.csv

Required Arguments:

  Input Data:
    --database_path       /Path/To/CSV.
Optional Arguments:
    --output_folder       Folder for output files (default $projectDir/results).

    """.stripIndent()
}


// Main workflow
workflow {


log.info """\
=============================================================================
         output_folder          : ${params.output_folder}
         database_path          : ${params.database_path}
         path_to_xml            : ${params.path_to_xml}
         num_points             : ${params.num_points}
         """
         .stripIndent()


    // Supported charge types : 'am1bcc', 'am1-mulliken', 'gasteiger', 'resp'
    // Show help message if the user specifies the --help flag at runtime
    // or if any required params are not provided
    if ( params.help || params.database_path == null){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    }

    if ( params.database_path ){
        // Define the input CSV file
        input_csv = file(params.database_path)

        // Create a channel with the CSV file
        csv_channel = channel.fromPath(input_csv)
        //vapor_systems = build_solvents(vapor_points.combine(path_to_xml))
        //path_to_database = Channel.fromPath( params.database_path )
        torch_model = train_model(csv_channel)
    } else {
        helpMessage()
    }
}