import java.util.*;

public class Graph {
    private int [][] adj;
    private int verticesNumber = 0;
    private int edgesNumber = 0;

    private boolean [] bipartiteColor;

    public Graph(int verticesNumber){
        if (verticesNumber >= 1) {
            this.verticesNumber = verticesNumber;
            adj = new int[verticesNumber][verticesNumber];
            bipartiteColor = new boolean[verticesNumber];
        }
        else
            try {
                throw new Exception("number of vertices can't be negative or zero");
            } catch (Exception e) {
                e.printStackTrace();
            }
    }

    public Graph(int [][] adj) {
        try {
            setAdj(adj);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void setAdj(int [][] adj) throws Exception {
        if (adj.length == adj[0].length){
            this.adj = adj;
            verticesNumber = adj.length;
            bipartiteColor = new boolean[verticesNumber];
            edgesNumber = getEdgesNumber();
        }
        else
            throw new Exception("row and column must be at the same size");
    }

    public int[][] getAdj() {
        return adj;
    }

    public int[][] getIncidentMatrix(){
        int edges = edgesNumber;
        int[][] incident = new int[verticesNumber][edges];
        int col = 0;
        if (isDirected()){
            for (int i = 0; i < verticesNumber; i++)
                for (int j = 0; j < verticesNumber; j++)
                    if (i == j) {
                        int l = col;
                        while (col < l + adj[i][j] / 2) {
                            incident[i][col] = 2;
                            col++;
                        }
                    }
                    else {
                        int m = col;
                        while (col < m + adj[i][j]) {
                            incident[i][col] = 1;
                            incident[j][col] = -1;
                            col++;
                        }
                    }
        }else
            for (int i = 0; i < verticesNumber - 1; i++)
                for (int j = i; j < verticesNumber; j++)
                    if (i == j) {
                        int l = col;
                        while (col < l + adj[i][j] / 2) {
                            incident[i][col] = 2;
                            col++;
                        }
                    }
                    else {
                        int m = col;
                        while (col < m + adj[i][j]) {
                            incident[i][col] = 1;
                            incident[j][col] = 1;
                            col++;
                        }
                    }
        return incident;
    }

    public void addEdge(int u, int v, boolean directed){
        if (!directed) addEdge(u, v);
        else adj[u][v]++;
    }

    public void addEdge(int u, int v){
        adj[u][v]++;
        adj[v][u]++;
    }


    public boolean isSimple() {
        for (int i = 0; i < verticesNumber; i++)
            for (int j = 0; j < verticesNumber; j++)
                if (adj[i][j] > 1) return false;
        return true;
    }

    public int getVerticesNumber(){
        return verticesNumber;
    }

    public int getEdgesNumber(){
        int sum = 0;
        if (isDirected()){
            for (int i = 0; i < verticesNumber; i++)
                for (int j = 0; j < verticesNumber; j++)
                    sum += adj[i][j];
        }
        else
            for (int i = 0; i < verticesNumber - 1; i++)
                for (int j = i; j < verticesNumber; j++)
                    sum += adj[i][j];

        return sum;
    }
    private int getEdgesNumber(int [][] a){
        int sum = 0;
        if (isDirected()){
            for (int i = 0; i < a.length; i++)
                for (int j = 0; j < a.length; j++)
                    sum += a[i][j];
        }
        else
            for (int i = 0; i < a.length - 1; i++)
                for (int j = i; j < a.length; j++)
                    sum += a[i][j];

        return sum;
    }

    public String determineType(){
        String s = "";
        if (isEmpty())
            return "empty";
        if (isSimple())
            s += " simple";
        if (isRegular())
            s += " regular";
        if (isComplete())
            s += " complete";
        if (isCycle())
            s += " cycle";
        if (isWheel())
            s += " wheel";
        if (isBipartite())
            s += " bipartite";
        if (isCompleteBipartite())
            s += " complete_bipartite";

        return s;
    }

    public boolean isDirected(){
        for (int i = 0; i < verticesNumber - 1; i++)
            for (int j = i + 1; j < verticesNumber; j++)
                if (adj[i][j] != adj[j][i]) return true;
        return false;
    }
    private boolean isCompleteBipartite() {
        if (isSimple() && isConnected() && isBipartite()) {
            int cFirst = 0;
            int cSec = 0;
            for (int i = 0; i < verticesNumber; i++)
                if (bipartiteColor[i]) cFirst++;
                else cSec++;

            return cFirst * cSec == edgesNumber;

//            for (int i = 0; i < verticesNumber; i++)
//                if (bipartiteColor[i]) {
//                    if (getVertexDeg(i) != cSec) return false;
//                } else if (getVertexDeg(i) != cFirst) return false;
//
//            return true;
        }
        return false;
    }

    public void printCycles(){
        ArrayList<ArrayList<Integer>> cycles = new ArrayList<>();
        for (int u = 0; u < verticesNumber; u++)
            for (int v = 0; v < verticesNumber; v++) {
                int[] color = new int[verticesNumber];
                int[] par = new int[verticesNumber];
                if (adj[u][v] > 0)
                    dfs_cycle(u, v, color, par, cycles);

                if (v >= u && adj[u][v] >= 2){
                    ArrayList<Integer> a = new ArrayList<>();
                    a.add(u);
                    a.add(v);
                    if (u != v)
                        a.add(u);
                    cycles.add(a);
                }
            }

        System.out.println(cycles.size());

        int count = 1;
        for (ArrayList<Integer> cycle: cycles) {
            System.out.println(count + "st:");
            System.out.println(cycle);
            count++;
        }
    }

    private void dfs_cycle(int u, int p, int [] color, int[] par, ArrayList<ArrayList<Integer>> cycles){
        if (color[u] == 2) return;

        if (color[u] == 1){
            int cur = p;
            ArrayList<Integer> c = new ArrayList<>();
            c.add(cur);
            while (cur != u){
                cur = par[cur];
                c.add(cur);
            }
            c.add(c.get(0));

            if(!checkCycleExists(cycles, c))
                cycles.add(c);

            return;
        }
        par[u] = p;

        color[u] = 1;

        for (int i = 0; i < verticesNumber; i++) {
            if (adj[u][i] > 0) {
                if (i != par[u]) {
                    dfs_cycle(i, u, color, par, cycles);
                }
            }
        }
        color[u] = 2;
    }

    private boolean checkCycleExists(ArrayList<ArrayList<Integer>> cycles, ArrayList<Integer> c) {
        for (ArrayList<Integer> cycle : cycles)
            if (Utils.compareLists(cycle, c))
                return true;
        return false;
    }

    private boolean containsLoop(){
        for (int i = 0; i < verticesNumber; i++)
            if (adj[i][i] != 0) return true;

        return false;
    }

    private boolean containsLoop(int [][] a){
        for (int i = 0; i < a.length; i++)
            if (a[i][i] != 0) return true;
        return false;
    }

    public boolean isBipartite() {
        if (isEmpty()) return false;
        if (containsLoop()) return false;
        LinkedList<Integer> q = new LinkedList<>();
        boolean [] visited = new boolean[verticesNumber];
        visited[0] = true;
        boolean [] color = new boolean[verticesNumber];
        boolean red = true;
        color[0] = red;
        q.add(0);
        while (q.size() != 0){
            int c = q.poll();
            red = !red;
            for (int i = 0; i < verticesNumber; i++){
                if (adj[c][i] > 0) {
                    if (!visited[i]) {
                        if (color[c] == red) return false;
                        color[i] = red;
                        visited[i] = true;
                        q.add(i);
                    } else
                    if (color[c] == color[i]) return false;
                }

            }
        }
        bipartiteColor = color;
        return true;
    }

    private boolean isWheel() {
        if (verticesNumber > 3 && isSimple() && isConnected()){
            if (verticesNumber == 4) return isComplete();
            boolean center = false;
            for (int i = 0; i < verticesNumber; i++)
                if (getVertexDeg(i) == verticesNumber - 1)
                    if (!center) center = true; else return false;
                else if (getVertexDeg(i) != 3) return false;

            return center;
        }

        return false;
    }

    private boolean isWheel(int [][] m) {
        //one element before last
        if (m[0][verticesNumber - 2] != 1) return false;
        if (m[verticesNumber - 2][0] != 1) return false;

        //upper diagonal
        for (int i = 0; i < verticesNumber - 1; i++)
            if (m[i][i+1] != 1) return false;

        //lower diagonal
        for (int i = 1; i < verticesNumber; i++)
            if (m[i][i-1] != 1) return false;

        //last row and column
        for (int i = 0; i < verticesNumber - 1; i++)
            if (m[i][verticesNumber - 1] != 1 && m[verticesNumber - 1][i] != 1)
                return false;

        return true;
    }

    private int[][] newIsomorphism(int first, int sec) {
        int [][] m = Utils.copyMatrix(adj);
        for (int i = 0; i < verticesNumber; i++){
            int temp = m[first][i];
            m[first][i] = m[sec][i];
            m[sec][i] = temp;
        }
        for (int i = 0; i < verticesNumber; i++){
            int temp = m[i][first];
            m[i][first] = m[i][sec];
            m[i][sec] = temp;
        }
        return m;
    }

    private int[][] newIsomorphism(int [][] a ,int first, int sec) {
        int [][] m = Utils.copyMatrix(a);
        for (int i = 0; i < verticesNumber; i++){
            int temp = m[first][i];
            m[first][i] = m[sec][i];
            m[sec][i] = temp;
        }
        for (int i = 0; i < verticesNumber; i++){
            int temp = m[i][first];
            m[i][first] = m[i][sec];
            m[i][sec] = temp;
        }
        return m;
    }

    public boolean isIsomorphic(Graph g){
        if (g.verticesNumber != verticesNumber || g.edgesNumber != edgesNumber) return false;
        if ((isSimple() && !g.isSimple()) || (g.isSimple() && !isSimple())) return false;
        if ((isConnected() && !g.isConnected()) || (g.isConnected() && !isConnected())) return false;
        if ((isCycle() && !g.isCycle()) || (g.isCycle() && !isCycle())) return false;

//        if ((isBipartite() && !g.isBipartite()) || (g.isBipartite() && !isBipartite())) return false;
        if ((isCompleteBipartite() && !g.isCompleteBipartite()) ||
                (g.isCompleteBipartite() && !isCompleteBipartite())) return false;
        if(!compareDegrees(g)) return false;
        if (g.isWheel() && isWheel()) return true;
        if (g.isRegular() && isRegular()) return true;
        return checkMatrices(g, 0);
    }

    private boolean checkMatrices(Graph g, int k) {
        if (g.compareGraphs(adj)) return true;

        for (int i = k; i < verticesNumber - 1; i++)
            for (int j = i + 1; j < verticesNumber; j++){
                int [][] m = newIsomorphism(g.getAdj(), i, j);
                if (checkMatrices(new Graph(m), k + 1)) return true;
            }

        return false;
    }

    private boolean compareDegrees(Graph g) {
        HashMap<Integer, Integer> h = new HashMap<>();
        HashMap<Integer, Integer> h2 = new HashMap<>();

        for (int i = 0; i < verticesNumber; i++){
            int d = getVertexDeg(i);
            int d2 = g.getVertexDeg(i);
            if (h.containsKey(d))
                h.replace(d, h.get(d) + 1);
            else
                h.put(d, 1);

            if (h2.containsKey(d2))
                h2.replace(d2, h2.get(d2) + 1);
            else
                h2.put(d2, 1);
        }
        for (Map.Entry e: h.entrySet())
            if (h2.containsKey(e.getKey())){
                if (h2.get(e.getKey()) != e.getValue())
                    return false;
            } else return false;

        return true;
    }

    public boolean compareGraphs(int[][] a) {
        for (int i = 0; i < verticesNumber; i++)
            for (int j = 0; j < verticesNumber; j++)
                if (a[i][j] != adj[i][j]) return false;
        return true;
    }

    private boolean isCycle() {
        if (verticesNumber == 1 && adj[0][0] == 2) return true;
        if (verticesNumber == 2 && adj[0][1] == 2 && adj[1][0] == 2) return true;

        if (isSimple() && isConnected()) {
            for (int i = 0; i < verticesNumber; i++)
                if (getVertexDeg(i) != 2)
                    return false;
            return true;
        }
        return false;
    }

    private boolean isComplete() {
        if (isSimple()) {
            for (int i = 0; i < verticesNumber; i++)
                if (getVertexDeg(i) != verticesNumber - 1)
                    return false;
            return true;
        }
        return false;
    }

    private boolean isRegular() {
        if (isSimple()) {
            int deg = getVertexDeg(0);
            for (int i = 1; i < verticesNumber; i++)
                if (deg != getVertexDeg(i))
                    return false;
            return true;
        }
        return false;
    }

    public void printComponents(){
        boolean [] visited = new boolean[verticesNumber];
        int count = 0;
        int x = 0;
        boolean all = true;
        while (all) {
            String s = "";
            LinkedList<Integer> q = new LinkedList<>();
            visited[x] = true;
            q.add(x);
            while (q.size() != 0) {
                int c = q.poll();
                s += " " + c;
                for (int i = 0; i < verticesNumber; i++) {
                    if (adj[c][i] > 0 && !visited[i]) {
                        visited[i] = true;
                        q.add(i);
                    }
                }
            }
            System.out.println(s);
            count++;
            for (int i = 0; i < verticesNumber; i++)
                if (!visited[i]) {
                    all = true;
                    x = i;
                    break;
                }
                else
                    all = false;
        }
        System.out.println(count + " components");
    }

    private boolean containsCycle2(int [][] a){
        if (containsLoop(a)) return true;
        return getEdgesNumber(a) >= a.length;
    }

    private boolean containsCycle(int [][] a){
        boolean [] visited = new boolean[a.length];
        for (int v = 0; v < a.length; v++){
            if (!visited[v])
                if (isCyclicUtil(v, visited, -1, a))
                    return true;
        }

        return false;
    }

    private boolean isCyclicUtil(int v, boolean[] visited, int p, int [][] a) {
        visited[v] = true;
        for (int u = 0; u < a.length; u++){
            if (a[v][u] > 0 )
                if (!visited[u]) {
                    if (isCyclicUtil(u, visited, v, a))
                        return true;
                } else if (u != p) return true;
        }
        return false;
    }

    public int spanningTreesNumber(){
        if (isComplete()) return (int) Math.pow(verticesNumber, verticesNumber - 2);
        if (isSimple()){
            int comEdges = verticesNumber * (verticesNumber - 1) / 2;
            if (edgesNumber == comEdges - 1) {
                int allSpanning = (int) Math.pow(verticesNumber, verticesNumber - 2);
                int spanningOfEdge = 2 * (int) Math.pow(verticesNumber, verticesNumber - 3);
                return allSpanning - spanningOfEdge;
            }
        }
        return getRecursiveSpanningTree(Utils.copyMatrix(adj));
    }

    private int getRecursiveSpanningTree(int [][] a) {
        if (isConnected(a)){
            if (containsLoop(a)) {
                int[][] t = removeLoops(Utils.copyMatrix(a));
                if (!containsCycle2(t))
                    return 1;
            }
            if (!containsCycle2(a))
                return 1;
        }
        else return 0;

        int [] v = chooseVertices(a);
        int [][] b = edgeReduction(Utils.copyMatrix(a), v[0], v[1]);
        int [][] c = edgeContraction(Utils.copyMatrix(a), v[0], v[1]);

        return getRecursiveSpanningTree(b) + getRecursiveSpanningTree(c);
    }

    private int[][] removeLoops(int[][] a) {
        for (int i = 0; i < a.length; i++)
            a[i][i] = 0;
        return a;
    }

    private int[] chooseVertices(int[][] a) {
        Random r = new Random();
        int f, l;
        do {
            f = r.nextInt(a.length);
            l = r.nextInt(a.length);
        }while (a[f][l] == 0 || f == l);
        return new int[]{f, l};
    }

    private int[][] edgeReduction(int[][] a, int f, int l) {
        a[f][l] = a[f][l] - 1;
        a[l][f] = a[l][f] - 1;
        return a;
    }

    private int[][] edgeContraction(int [][] a, int f, int l){
        int [][] b = new int[a.length - 1][a.length - 1];
        edgeReduction(a, f, l);
        for (int i = 0; i < a.length; i++)
            a[f][i] += a[l][i];
        for (int i = 0; i < a.length; i++)
            a[i][f] += a[i][l];

        int row = 0;
        for (int i = 0; i < a.length; i++) {
            if (i != l) {
                int col = 0;
                for (int j = 0; j < a.length; j++)
                    if (j != l) {
                        b[row][col] = a[i][j];
                        col++;
                    }
                row++;
            }
        }

        return b;
    }

    public int getRadius(){
        int [][] c = allPairsShortest();

        int [] max = new int[verticesNumber];
        for (int i =0; i < verticesNumber; i++) {
            max[i] = c[i][0];
            for (int j = 1; j < verticesNumber; j++)
                if (c[i][j] > max[i]) max[i] = c[i][j];
        }
        int min = max[0];
        for (int i = 1; i < verticesNumber; i++)
            if (max[i] < min) min = max[i];

        return min;
    }

    public int getCenter(){
        int [][] c = allPairsShortest();

        int [] max = new int[verticesNumber];
        for (int i =0; i < verticesNumber; i++) {
            max[i] = c[i][0];
            for (int j = 1; j < verticesNumber; j++)
                if (c[i][j] > max[i]) max[i] = c[i][j];
        }

        int min = max[0];
        int v = 0;
        for (int i = 1; i < verticesNumber; i++)
            if (max[i] < min) {
                min = max[i];
                v = i;
            }
        return v;
    }

    private int[][] allPairsShortest() {
        int [][] c = Utils.copyMatrix(adj);
        for (int i = 0; i < verticesNumber; i++)
            for (int j = 0; j < verticesNumber; j++)
                if (i != j && c[i][j] == 0) c[i][j] = 99;

        for (int k = 0; k < verticesNumber; k++)
            for (int i = 0; i < verticesNumber; i++)
                for (int j = 0; j < verticesNumber; j++)
                    if (c[i][k] + c[k][j] < c[i][j])
                        c[i][j] = c[i][k] + c[k][j];

        return c;
    }

    private boolean isConnected(int [][] a){
        LinkedList<Integer> q = new LinkedList<>();
        boolean [] visited = new boolean[a.length];
        visited[0] = true;
        q.add(0);
        while (q.size() != 0){
            int c = q.poll();

            for (int i = 0; i < a.length; i++){
                if (a[c][i] > 0 && !visited[i]){
                    visited[i] = true;
                    q.add(i);
                }
            }
        }
        for (int i = 0; i < a.length; i++)
            if (!visited[i]) return false;
        return true;
    }

    public boolean isConnected(){
        LinkedList<Integer> q = new LinkedList<>();
        boolean [] visited = new boolean[verticesNumber];
        visited[0] = true;
        q.add(0);
        while (q.size() != 0){
            int c = q.poll();

            for (int i = 0; i < verticesNumber; i++){
                if (adj[c][i] > 0 && !visited[i]){
                    visited[i] = true;
                    q.add(i);
                }
            }
        }
        for (int i = 0; i < verticesNumber; i++)
            if (!visited[i]) return false;
        return true;
    }

    public int getVertexDeg(int vertex) {
        int sum = 0;
        for (int i = 0 ; i < verticesNumber; i++)
            sum += adj[i][vertex];
        return sum;
    }

    private boolean isEmpty() {
        for (int i = 0; i < verticesNumber; i++)
            for (int j = 0; j < verticesNumber; j++)
                if (adj[i][j] != 0) return false;
        return true;
    }

}