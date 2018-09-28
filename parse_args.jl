using ArgParse: ArgParseSettings, @add_arg_table, parse_args

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
            default = 2321499335
        "--chromosome-sizes"
            help = "Sizes of the chromosomes. Only needed when no genome is given."
            dest_name = "chromosome_sizes"
        "--outfile"
            help = "File to write results to."
            range_tester = create_path_or_fail
            required = true
    end

    return parse_args(s)
end


function main()

    parsed_args = parse_commandline()

    println("Parsed args:")
    println(parsed_args)

    for (arg, val) in parsed_args
        if val != nothing
            println("  $arg => $val")
        else
            println("  Missing: $arg")
        end
    end
end

# main()
