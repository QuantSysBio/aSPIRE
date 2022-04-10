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
		qValCutoff = config["qVal"]
	conda:
		"dependencies.yaml"
	script:
		"createInput.R"


rule wait_for_skyline:
	input:
		ssl = "results/{protein_name}/{protein_name}.ssl",
		fasta = "results/{protein_name}/{protein_name}.fasta"
	output:
		skyline = "results/{protein_name}/MS1_HPR.csv"
	params:
		protein_name = config["protein_name"]
	conda:
		"dependencies.yaml"
	script:
		"wait.R"


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
		plot_raw = "results/{protein_name}/plots/rawIntensities_noSynErr.pdf",
		plot_err = "results/{protein_name}/plots/rawIntensities_SynErr.pdf",
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
		plot_norm = "results/{protein_name}/plots/normIntensities_allReps.pdf",
		plot_kinetics = "results/{protein_name}/plots/finalKinetics.pdf",
		final_kinetics = "results/{protein_name}/finalKinetics.csv"
	params:
		protein_name = config["protein_name"]
	conda:
		"dependencies.yaml"
	script:
		"postprocessing.R"




