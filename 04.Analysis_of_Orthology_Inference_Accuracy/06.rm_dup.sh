awk -F, 'BEGIN{OFS=","}
NR==1{
    # Record header
    for(i=1;i<=NF;i++) hdr[i]=$i
    next
}
{
    key=$1
    if(!(key in seen)){ order[++n]=key; seen[key]=1 }
    for(i=2;i<=NF;i++){
        # Logical OR across rows = take the maximum value
        if((key FS i) in val){
            if($i>val[key FS i]) val[key FS i]=$i+0
        } else {
            val[key FS i]=$i+0
        }
    }
}
END{
    # Print header
    for(i=1;i<=length(hdr);i++){
        printf "%s%s", hdr[i], (i<length(hdr)?OFS:ORS)
    }
    # Output per species in original order of appearance
    for(k=1;k<=n;k++){
        key=order[k]
        printf "%s", key
        for(i=2;i<=length(hdr);i++){
            printf "%s%d", OFS, val[key FS i]+0
        }
        printf "\n"
    }
}' 37-ACEK-202_0-1.transposed.csv > 37-ACEK-202_final.csv

# Identify numeric columns (from column 2 onward) whose sum equals 0
# and write their names to zero_columns.txt
awk -F, 'NR==1{
    for(i=1;i<=NF;i++) hdr[i]=$i
    next
}
{
    for(i=2;i<=NF;i++){
        sum[i]+=$i+0
    }
}
END{
    cnt=0
    for(i=2;i in sum;i++){
        if(sum[i]==0){
            zero[++cnt]=hdr[i]
        }
    }
    # Print to terminal
    print "columns_with_sum_zero_count=" cnt
    if(cnt>0){
        print "columns_with_sum_zero:"
        for(i=1;i<=cnt;i++) print zero[i]
    }
    # Also write to file
    out="zero_columns.txt"
    print "columns_with_sum_zero_count=" cnt > out
    if(cnt>0){
        print "columns_with_sum_zero:" >> out
        for(i=1;i<=cnt;i++) print zero[i] >> out
    }
}' 37-ACEK-202_final.csv