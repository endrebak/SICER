using ArgParse: ArgParseSettings, @add_arg_table, parse_args

using BioAlignments


function find_closest_readlength(file)

    df = nothing
    try
        df = file_to_df(file, 9999)
    catch # in case file less than 10k lines
        df = file_to_df(file)
    end

    lengths = df[3] - df[2]
    estimated_readlength = Statistics.median(df[3] - df[2])

    readlengths = [36, 50, 75, 100]
    differences = [abs(r - estimated_readlength) for r in readlengths]
    min_difference = minimum(differences)
    index_of_min_difference = [i
                               for (i, d) in enumerate(differences)
                               if d == min_difference][1]

    readlengths[index_of_min_difference]

end


function read_egf(genome, closest_readlength)
    datapath = joinpath(@__DIR__, "genomes/")
    f = joinpath(datapath, "effective_genome_sizes", lowercase(genome) * string("_") * string(closest_readlength) * ".txt")
    df = CSV.read(f, delim="\t", header=false)
    parse(Float64, split(df[4, 1], ":  ")[2])
end


function find_egf(file, genome)
    closest_read_length = find_closest_readlength(file)
    read_egf(genome, closest_read_length)
end

function chromosome_sizes(genome)
    # datapath = joinpath(@__DIR__, "genomes/")
    datapath = joinpath(@__DIR__, "genomes/")
    # println(datapath)
    f = joinpath(datapath, "chromsizes", lowercase(genome) * ".chromsizes")
    # println(f)
    CSV.read(f, delim="\t", header=false)
end



function create_path_or_fail(x)

    # will raise SystemError if trying to create a file which
    # we are not allowed to create. e.g. ("/usr/hello.txt")
    result = mkpath(dirname(x))
    open(x, "w+")

    return true # isdir(result) || result == nothing
end

function parse_commandline()
    s = ArgParseSettings("SICER.jl",
                         version = "Version 0.0.1")

    @add_arg_table s begin
        "--chip"
            help = "One or more ChIP files to analyze"
            nargs = '+'
            range_tester = isfile # File not found
        "--input", "-o"
            help = "Zero or more Input files to use as background"
            nargs = '*'
            range_tester = isfile # File not found
        "--keep-duplicates"
            help = "Keep all duplicate reads. Default is to remove all but first duplicate read."
            action = :store_true
            dest_name = "keep_duplicates"
        "--fragment-size"
            help = "Size of fragments."
            default = 150
            dest_name = "fragment_size"
            arg_type=Int64
        "--bin-size"
            help = "Size of bins to count reads in."
            arg_type = Int64
            default = 200
            dest_name = "bin_size"
        "--gaps_allowed"
            help = "Number of gaps allowed."
            arg_type = Int64
            default = 3
            dest_name = "gaps_allowed"
        "--e-value"
            help = "E-value for identification of significant islands. Only used when no --input files are given."
            default = 1000
            arg_type=Int64
            dest_name = "e_value"
        "--genome"
            help = "The genome to use"
            default = "hg19"
            dest_name = "genome"
        "--effective-genome-fraction"
            help = "The effective genome fraction to use. Only needed when no genome is given."
            dest_name = "effective_genome_fraction"
            default = 0.80
        "--chromosome-sizes"
            help = "Sizes of the chromosomes. Only needed when no genome is given."
            dest_name = "chromosome_sizes"
        "--outfile"
            help = "File to write results to."
            range_tester = create_path_or_fail
            required = true
    end

    args = parse_args(s)

    if args["genome"] != nothing
        genome = lowercase(args["genome"])
        println("find egf")
        args["effective_genome_fraction"] = find_egf(args["chip"][1], genome)
        println("find chromsizes")
        args["chromosome_sizes"] = chromosome_sizes(genome)
        args["effective_genome_size"] = sum(args["chromosome_sizes"][2]) * args["effective_genome_fraction"]
    else
        args["chromosome_sizes"] = CSV.read(f, delim="\t", header=false)
        args["effective_genome_size"] = sum(args["chromosome_sizes"][2]) * args["effective_genome_fraction"]
    end

    args
end


function main()

    parsed_args = parse_commandline()

    for (arg, val) in parsed_args
        if val != nothing
            println("  $arg => $val")
        else
            println("  Missing: $arg")
        end
    end
end
