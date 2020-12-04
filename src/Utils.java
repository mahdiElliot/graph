import java.util.ArrayList;

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
}
