# utilities.jl — Low-level parsing helpers for NASTRAN BDF format

function safe_get(arr::Vector{Any}, idx::Int, default_val=nothing)
    if idx > length(arr); return default_val; end
    return arr[idx]
end

function parse_nastran_number(field::Any, default_val=nothing)
    if isa(field, Number); return field; end
    if isnothing(field) || strip(string(field)) == ""; return default_val; end

    field_str = String(strip(string(field)))
    clean_field = replace(field_str, r"([\d.])([+-])(\d)" => s"\1e\2\3")
    clean_field = replace(clean_field, "ee" => "e")

    try
        val = parse(Float64, clean_field)
        if abs(val) >= 0.5 && abs(val - round(val)) < 1e-8; return Int(round(val)); end
        return val
    catch; return default_val; end
end

function to_id(val)
    if isa(val, Integer); return val; end
    if isa(val, AbstractFloat); return Int(round(val)); end
    return 0
end

function expand_nastran_list(raw_fields)
    result = Int[]
    i = 1
    while i <= length(raw_fields)
        val = raw_fields[i]
        if isa(val, AbstractString) && uppercase(strip(val)) == "THRU"
            if isempty(result) || i == length(raw_fields)
                i += 1; continue
            end
            start_id = result[end]
            end_id = to_id(parse_nastran_number(raw_fields[i+1]))
            if end_id > start_id; append!(result, (start_id+1):end_id); end
            i += 2
        else
            parsed = to_id(parse_nastran_number(val, 0))
            if parsed > 0; push!(result, parsed); end
            i += 1
        end
    end
    return result
end

function get_nastran_card_name(line::AbstractString)
    if occursin(",", line)
        parts = split(line, ",")
        return uppercase(strip(parts[1]))
    end
    return uppercase(strip(line[1:min(8, length(line))]))
end

function get_nastran_fields_from_line(line::AbstractString)
    fields = []
    if occursin(",", line)
        parts = split(line, ",")
        if length(parts) > 1
            for p in parts[2:end]; push!(fields, strip(string(p))); end
        end
        # Pad free-field lines to 8 fields to maintain NASTRAN card field alignment
        while length(fields) < 8; push!(fields, ""); end
    else
        padded = rpad(line, 80, ' ')
        for i in 0:7
             s = 9 + i*8; e = s + 7
             if s > length(padded); break; end
             push!(fields, strip(padded[s:min(e, end)]))
        end
    end
    return fields
end
