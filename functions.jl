"""
    string_sites!(df)

Convert the `site` column of `df` to strings, so site ids read by CSV stay strings.
Returns `df`.
"""
function string_sites!(df)
    df.site = string.(df.site)
    return df
end
