package utils;
import java.util.*;

public class ArrayList2<Item> extends ArrayList<Item> {
    public ArrayList2(Collection<Item> arr){
        super(arr);
    }

    public ArrayList<ArrayList<Item>> subsets(int filterSubsetSize){
        ArrayList<ArrayList<Item>> arr = new ArrayList<>();
        int n = this.size();
        for (int i = 0; i < (1 << n); i++){
            ArrayList<Item> s = new ArrayList<>();
            for (int j = 0; j < n; j++)
                if ((i & (1 << j)) > 0)
                    s.add(this.get(j));

            if (s.size() >= filterSubsetSize) arr.add(s);
        }
        return arr;
    }
}
