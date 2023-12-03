#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Using recursion
nextflow.preview.recursion=true

// All of the default parameters are being set in `nextflow.config`

// Import sub-workflows
include { build_solvents } from './modules/system_builder'
include { train_model } from './modules/model_builder'
include { initialize_scikit_optimize_model } from './modules/scikit_optimize'
include { calibrate } from './modules/scikit_optimize'



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
         batch_size             : ${params.batch_size}
         density_lb             : ${params.density_lb}
         density_ub             : ${params.density_ub}
         alpha_lb               : ${params.alpha_lb}
         alpha_ub               : ${params.alpha_ub}
         torch_model            : ${params.torch_model}
         torch_scalers          : ${params.torch_scalers}


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
        /**
        train_model(csv_channel)
        params.torch_model = train_model.out.torch_model
        params.torch_scalers = train_model.out.torch_scalers
        train_model.out.torch_model.view()
        println("torch_model")
        println("${params.torch_model}")
        println("torch_scalers")
        println("${params.torch_scalers}")
        */
        skopt_model = initialize_scikit_optimize_model()
        calibrate.recurse(skopt_model).times(3)
    } else {
        helpMessage()
    }
}