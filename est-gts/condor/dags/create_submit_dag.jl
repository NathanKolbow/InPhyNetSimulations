cd("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/condor/dags/")


ntaxa = parse(Int, ARGS[1])
rep = parse(Int, ARGS[2])
ils = ARGS[3]
ngt = parse(Int, ARGS[4])
m = parse(Int, ARGS[5])
resubmit = false
if length(ARGS) == 6 && ARGS[6] == "true"
    resubmit = true
end


# 1. create directory for dag
dag_dir = "$(pwd())/n$(ntaxa)"
if !isdir(dag_dir) mkdir(dag_dir) end
dag_dir *= "/$(ils)-$(ngt)gt-$(m)m"
if !isdir(dag_dir) mkdir(dag_dir) end
dag_dir *= "/rep$(rep)"
if !isdir(dag_dir) mkdir(dag_dir) end
if !isdir(joinpath(dag_dir, "err")) mkdir(joinpath(dag_dir, "err")) end
if !isdir(joinpath(dag_dir, "out")) mkdir(joinpath(dag_dir, "out")) end


# 2. Copy template
dag_file = "$(dag_dir)/dag.dag"
if !isfile(dag_file)
    cp("dag_template.dag", dag_file)
    lines = readlines(dag_file)

    open(dag_file, "w+") do f
        for line in lines
            new_line = replace(line, "[[NTAXA]]" => "$(ntaxa)",
                                    "[[REP]]" => "$(rep)",
                                    "[[ILS]]" => ils,
                                    "[[NGT]]" => "$(ngt)",
                                    "[[M]]" => "$(m)",
                                    "[[SUB_DIR]]" => "/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/condor/dags/submit_files")
            write(f, "$(new_line)\n")
        end
    end
    printstyled("[CREATED] ", color=:green)
    println("$(dag_file)")


    printstyled("[SUBMITTING]", color=:yellow)
    cd(dag_dir)
    run(`condor_submit_dag $(dag_file)`)
else
    printstyled("[ALREADY EXISTS] ", color=:cyan)
    if !resubmit
        println("DAG file already exists. Not resubmitting.")
    elseif (ntaxa, rep, ils, ngt, m) != (500, 2, "high", 1000, 10) && (ntaxa, rep, ils, ngt, m) != (500, 2, "high", 1000, 20) && (ntaxa, rep, ils, ngt, m) != (500, 3, "high", 1000, 10) && (ntaxa, rep, ils, ngt, m) != (500, 4, "low", 1000, 20) && (ntaxa, rep, ils, ngt, m) != (500, 2, "low", 1000, 20) && (ntaxa, rep, ils, ngt, m) != (500, 5, "low", 100, 10)
        println("Resubmitting DAG.")
        run(`condor_submit_dag -f $(dag_file)`)
    end
end
