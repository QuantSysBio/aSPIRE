wildcard_constraints:
	protein_name = config["protein_name"]


rule create_input:
	input:
		sample_list = "data/sample_list.csv"
	output:
		assignments_csv = "results/{protein_name}/ASSIGNMENTS.csv",
		assignments_bin = "results/{protein_name}/ASSIGNMENTS.RData",
		ssl = "results/{protein_name}/{protein_name}.ssl",
		fasta = "results/{protein_name}/{protein_name}.fasta"
	params:
		protein_name = config["protein_name"],
		spAngleCutoff = config["spAngle"],
		qValCutoff = config["qVal"],
		rtCutoff = config["RT"]
	conda:
		"dependencies.yaml"
	script:
		"createInput.R"


rule run_skyline:
	input:
		ssl = "results/{protein_name}/{protein_name}.ssl",
		fasta = "results/{protein_name}/{protein_name}.fasta"
	output:
		skyline = "results/{protein_name}/MS1_HPR.csv",
		tics = "results/{protein_name}/{protein_name}_TICs.tsv"
	params:
		protein_name = config["protein_name"],
		raw_file_loc = config["raw_file_loc"],
		automatedSkyline = config["automatedSkyline"]
	conda:
		"dependencies.yaml"
	script:
		"Skyline.R"


rule get_quantities:
	input:
		sample_list = "data/sample_list.csv",
		skyline = "results/{protein_name}/MS1_HPR.csv",
		assignments_bin = "results/{protein_name}/ASSIGNMENTS.RData"
	output:
		quantities_raw = "results/{protein_name}/QUANTITIES_raw.RData"
	params:
		protein_name = config["protein_name"]
	conda:
		"dependencies.yaml"
	script:
		"getQuant.R"


rule map_and_filter:
	input:
		sample_list = "data/sample_list.csv",
		assignments_bin = "results/{protein_name}/ASSIGNMENTS.RData",
		quantities_raw = "results/{protein_name}/QUANTITIES_raw.RData"
	output:
		quantities_filtered = "results/{protein_name}/QUANTITIES_filtered.RData"
	params:
		protein_name = config["protein_name"]
	conda:
		"dependencies.yaml"
	script:
		"mapRemSynErr.R"


rule postprocessing:
	input:
		quantities_filtered = "results/{protein_name}/QUANTITIES_filtered.RData",
		assignments_bin = "results/{protein_name}/ASSIGNMENTS.RData"
	output:
		plot_raw = "results/{protein_name}/plots/{protein_name}_rawIntensities.pdf",
		plot_kinetics = "results/{protein_name}/plots/{protein_name}_finalKinetics.pdf",
		final_kinetics = "results/{protein_name}/{protein_name}_finalKinetics.csv"
	params:
		protein_name = config["protein_name"]
	conda:
		"dependencies.yaml"
	script:
		"postprocessing.R"


rule plot_tics:
	input:
		sample_list = "data/sample_list.csv",
		tics = "results/{protein_name}/{protein_name}_TICs.tsv"
	output:
		chromatogram = "results/{protein_name}/plots/{protein_name}_TICs.done"
	params:
		protein_name = config["protein_name"]
	conda:
		"dependencies.yaml"
	script:
		"plotTICs.R"


rule plot_coverage:
	input:
		final_kinetics = "results/{protein_name}/{protein_name}_finalKinetics.csv"
	output:
		coveragevalues = "results/{protein_name}/plots/{protein_name}_coveragevalues.RData",
		residuemap = "results/{protein_name}/plots/{protein_name}_residuemap.gif"
	params:
		protein_name = config["protein_name"]
	conda:
		"dependencies.yaml"
	script:
		"coverage_maps.R"

