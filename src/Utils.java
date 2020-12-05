import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

public final class Utils {
    public static void printMatrix(int[][] m) {
        String s = "";
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[i].length; j++)
                s += m[i][j] + " ";
            s += "\n";
        }
        System.out.println(s);
    }

    static int[][] copyMatrix(int[][] m) {
        int[][] b = new int[m.length][m.length];
        for (int i = 0; i < m.length; i++)
            for (int j = 0; j < m[i].length; j++)
                b[i][j] = m[i][j];

            return b;
    }

    static void printArray(int [] arr){
        String s = arr[0] + "";
        for (int i = 1; i < arr.length; i++)
           s += " " + arr[i];

        System.out.println(s);
    }

    static boolean compareLists(ArrayList<Integer> l1, ArrayList<Integer> l2){
        if (l1.size() != l2.size()) return false;

        for (int l: l1)
            if (!l2.contains(l)) return false;

        return true;
    }

    static int fact(int n){
        if (n == 1 || n == 0)
            return 1;
        return n * fact(n-1);
    }

    static ArrayList<ArrayList<Integer>> subsets(int n, int filterSubsetSize){
        ArrayList<ArrayList<Integer>> arr = new ArrayList<>();
        for (int i = 1; i < (1<<n); i++) {
            ArrayList<Integer> a = new ArrayList<>();
            // Print current subset
            for (int j = 0; j < n; j++)
                // (1<<j) is a number with jth bit 1
                // so when we 'and' them with the
                // subset number we get which numbers
                // are present in the subset and which
                // are not
                if ((i & (1 << j)) > 0)
                    a.add(j);

            if (a.size() > filterSubsetSize)  arr.add(a);
        }
        return arr;
    }

    static ArrayList<ArrayList<Edge>> subsets(ArrayList<Edge> set, int filterSubsetSize){
        ArrayList<ArrayList<Edge>> arr = new ArrayList<>();
        int n = set.size();
        for (int i = 0; i < (1 << n); i++){
            ArrayList<Edge> s = new ArrayList<>();
            for (int j = 0; j < n; j++)
                if ((i & (1 << j)) > 0)
                    s.add(set.get(j));

            if (s.size() > filterSubsetSize) arr.add(s);
        }

        return arr;
    }
}
