import java.util.ArrayList;
import java.util.*;

public class Test {
    public static void main(String[] args) {
//        int [][] adj = {{0, 1, 1, 0, 0, 1},
//                        {1, 0, 1, 1, 0, 0},
//                        {1, 1, 0, 1, 1, 1},
//                        {0, 1, 1, 0, 1, 0},
//                        {0, 0, 1, 1, 0, 1},
//                        {1, 0, 1, 0, 1, 0}};

//        int [][] adj = {
//                {0, 1, 0, 0},
//                {1, 0, 0, 1},
//                {0, 0, 0, 0},
//                {0, 1, 0, 0}
//        };

//        int [][] adj = {
//                {0, 1, 0, 1},
//                {1, 0, 1, 0},
//                {0, 1, 0, 1},
//                {1, 0, 1, 0}
//        };
//        int [][] adj = {
//                {0, 1, 1},
//                {1, 0, 1},
//                {1, 0, 0}
//        };

//        int [][] adj = {
//                {0, 1, 0, 0},
//                {0, 0, 1, 0},
//                {0, 0, 0, 1},
//                {1, 0, 0, 0}
//        };

//        int[][] adj = {
//                {0, 1, 1, 0, 0, 0},
//                {1, 0, 1, 0, 0, 0},
//                {1, 1, 0, 0, 0, 0},
//                {0, 0, 0, 0, 1, 1},
//                {0, 0, 0, 1, 0, 1},
//                {0, 0, 0, 1, 1, 0}
//        };

//        int [][] adj = {
//                {0, 1, 0, 0, 1},
//                {1, 0, 1, 1, 0},
//                {0, 1, 0, 1, 0},
//                {0, 1, 1, 0, 1},
//                {1, 0, 0, 1, 0}
//        };

//        int [][] adj = {
//                {0, 1, 0, 0, 0, 0, 0, 0, 0, 1},
//                {1, 0, 1, 1, 1, 0, 0, 0, 0, 1},
//                {0, 1, 0, 1, 0, 0, 0, 0, 0, 0},
//                {0, 1, 1, 0, 0, 0, 0, 0, 0, 0},
//                {0, 1, 0, 0, 0, 1, 0, 0, 1, 0},
//                {0, 0, 0, 0, 1, 0, 1, 1, 0, 0},
//                {0, 0, 0, 0, 0, 1, 0, 1, 0, 0},
//                {0, 0, 0, 0, 0, 1, 1, 0, 0, 0},
//                {0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
//                {1, 1, 0, 0, 0, 0, 0, 0, 0, 0}
//        };

//        int [][] adj = {
//                {0, 1, 0, 0},
//                {1, 0, 1, 0},
//                {0, 1, 0, 1},
//                {0, 0, 1, 0}
//        };

//        int [][] adj = {
//                {0, 0, 0, 1},
//                {0, 0, 1, 0},
//                {0, 1, 0, 1},
//                {1, 0, 1, 0}
//        };

//        int [][] adj = {
//                {0, 1, 2},
//                {1, 0, 0},
//                {2, 0, 0}
//        };

//        int [][] adj = {
//                {0, 1, 0, 1, 1},
//                {1, 0, 1, 0, 1},
//                {0, 1, 0, 1, 1},
//                {1, 0, 1, 0, 1},
//                {1, 1, 1, 1, 0}
//        };

//        int [][] adj = {
//                {0, 2, 0, 3},
//                {2, 0, 2, 1},
//                {0, 2, 0, 4},
//                {3, 1, 4, 0}
//        };

//        int [][] adj = {
//                {0, 1, 1, 1, 1},
//                {1, 0, 1, 1, 1},
//                {1, 1, 0, 1, 1},
//                {1, 1, 1, 0, 1},
//                {1, 1, 1, 1, 0}
//        };

//        int [][] adj = {
//                {0, 1, 1},
//                {1, 0, 1},
//                {1, 1, 0}
//        };

        int [][] adj = {
                {0, 1, 0, 0, 1},
                {1, 0, 1, 0, 0},
                {0, 1, 0, 1, 1},
                {0, 0, 1, 0, 1},
                {1, 0, 1, 1, 0}
        };

//        int [][] adj = {
//                /*a1*/{0, 1, 1, 0, 0, 0},
//                /*a2*/{1, 0, 1, 1, 0, 0},
//                /*a3*/{1, 1, 0, 1, 0, 0},
//                /*a4*/{0, 1, 1, 0, 1, 0},
//                /*a5*/{0, 0, 0, 1, 0, 0},
//                /*a6*/{0, 0, 0, 0, 0, 0},
//        };

//        int [][] adj = {
//                {0, 1, 0, 0, 0, 1, 0},
//                {1, 0, 1, 0, 0, 1, 1},
//                {0, 1, 0, 1, 1, 0, 1},
//                {0, 0, 1, 0, 1, 0, 0},
//                {0, 0, 1, 1, 0, 1, 1},
//                {1, 1, 0, 0, 1, 0, 1},
//                {0, 1, 1, 0, 1, 1, 0}
//        };

//        int [][] adj = {
//                {0, 1, 0, 1},
//                {1, 0, 1, 0},
//                {0, 1, 0, 1},
//                {1, 0, 1, 0}
//        };

//        int [][] adj = {
//                {0, 1},
//                {1, 0}
//        };
//        int [][] adj2 = {
//                {0, 0, 0, 1},
//                {0, 0, 1, 1},
//                {0, 1, 0, 0},
//                {1, 1, 0, 0}
//        };

        int [][] adj2 = {
                {0, 1, 1, 0},
                {1, 0, 0, 1},
                {1, 0, 0, 1},
                {0, 1, 1, 0}
        };

        Graph g = new Graph(adj);
        Graph g2 = new Graph(adj2);
//        System.out.println(g.isDirected());
//        System.out.println(g.maxIndependentSet());
//        System.out.println(g.minEdgeCover());
//        System.out.println(g.maxMatchingNumber());
//        utils.Utils.printMatrix(g.getIncidentMatrix());
//        System.out.println(g.determineType());
//        System.out.println(g.isConnected());
//        System.out.println(g.isDirected());
//        g.getComponents();
//        System.out.println(g.getCenter());
//        System.out.println(g.spanningTreesNumber());
//        System.out.println(g.containsCycle(adj));
//        System.out.println(g.isIsomorphic(g2));
//        g.printCycles();

//        System.out.println(g.eulerianWalk());

//        boolean [] bipartiteColor = {true, false, true, false};
//        System.out.println(g.checkHall(true, bipartiteColor));
//        System.out.println(g.isBipartite2());
        System.out.println(g.numberOfComponents());
//        System.out.println(g.isIsomorphic2(g2));
    }
}

