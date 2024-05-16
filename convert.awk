# AWK script to read a CSV file and convert float values to 2 digits after the dot

BEGIN {
    FS = ","; # Set field separator as comma for CSV
    OFS = ","; # Set output field separator as comma
}

{
    for (i = 1; i <= NF; i++) { # Loop through all fields in a record
        if ($i ~ /^[+-]?[0-9]*[.][0-9]+$/) { # Check if field is a float
            $i = sprintf("%.1f", $i); # Convert float to 2 digits after the dot
        }
    }
    print; # Print the modified record
}