for pop in temp trop
do
    echo -e "Working on ${pop}\n"
    head -1 window_nd/window_10:100010913-100013909_${pop}.txt > window_nd_${pop}.txt
    find ./window_nd/ -name "*${pop}.txt" | xargs tail -qn +2 >> window_nd_${pop}.txt
done

    
    
