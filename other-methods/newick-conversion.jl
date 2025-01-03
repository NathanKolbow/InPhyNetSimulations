function remove_nested_markers(text::String)::String
    # Sub 1
    # Define the regex pattern for markers
    pattern1 = r"<\d+\|(.*)"

    # Iteratively remove markers until no more matches
    while occursin(pattern1, text)
        text = replace(text, pattern1 => s"\1")
    end

    # Sub 2
    # Define the regex pattern for markers
    pattern2 = r"(.*)\|\d+:\d+\.\d+>"

    # Iteratively remove markers until no more matches
    while occursin(pattern2, text)
        text = replace(text, pattern2 => s"\1")
    end

    # Sub 3
    # Define the regex pattern for markers
    pattern3 = r"(.*)\|\d+>"

    # Iteratively remove markers until no more matches
    while occursin(pattern3, text)
        text = replace(text, pattern3 => s"\1")
    end

    return text
end