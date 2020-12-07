import jdk.jshell.execution.Util;

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
            setEdgesNumber();
        }
        else
            throw new Exception("row and column must be at the same size");
    }

    private void setEdgesNumber() {
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

        this.edgesNumber = sum;
    }

    public int getEdgesNumber(){
        return edgesNumber;
    }

    public int[][] getAdj() {
        return adj;
    }

    public int[][] getIncidentMatrix(){
        int[][] incident = new int[verticesNumber][edgesNumber];
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
        else {
            adj[u][v]++;
            edgesNumber++;
        }
    }

    public void addEdge(int u, int v){
        adj[u][v]++;
        adj[v][u]++;
        edgesNumber++;
    }

    public boolean removeEdge(int u, int v, boolean directed) {
        if (!directed) return removeEdge(u, v);
        else if (adj[u][v] > 0) {
            adj[u][v]--;
            edgesNumber--;
            return true;
        }
        return false;
    }

    private boolean removeEdge(int u, int v){
        if (adj[u][v] > 0) {
            adj[u][v]--;
            adj[v][u]--;
            edgesNumber--;
            return true;
        }
        return false;
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
        if (isPath())
            s += " path";
        return s;
    }

    private boolean isPath(){
        if (verticesNumber == 1 && edgesNumber == 0) return true;
        if (isDirected()){
            if (verticesNumber == 2 && edgesNumber == 1 &&
                    ((adj[0][1] == 1 && adj[1][0] == 0) || (adj[0][1] == 0 && adj[1][0] == 1))) return true;
            if (isSimple() && !containsCycle(adj)){
                int number = 0;
                boolean end = false;
                for (int i = 0; i < verticesNumber; i++){
                    int d = getVertexDegree(i);
                    if (d == 1) number++;
                    else if (d == 0) end = true;
                    else return false;
                }
                return end && number == verticesNumber - 1;
            }
        }
        else {
            if (verticesNumber == 2 && edgesNumber == 1 && adj[0][1] == 1 && adj[1][0] == 1) return true;
            if (isSimple() && isConnected() && !containsCycle(adj)){
                short numberOfEndpoints = 0;
                for (int i = 0; i < verticesNumber; i++){
                    int d = getVertexDegree(i);
                    if (d == 1)
                        numberOfEndpoints++;
                    else if (d > 2)
                        return false;
                    if (numberOfEndpoints > 2) return false;
                }
                return true;
            }
        }
        return false;
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
        }
        return false;
    }

    private boolean isCompleteRegularBipartite(){
        if (isSimple() && isConnected() && isBipartite()) {
            int cFirst = 0;
            int cSec = 0;
            for (int i = 0; i < verticesNumber; i++)
                if (bipartiteColor[i]) cFirst++;
                else cSec++;

            return cFirst * cSec == edgesNumber && cFirst == cSec;
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
                if (getVertexDegree(i) == verticesNumber - 1)
                    if (!center) center = true; else return false;
                else if (getVertexDegree(i) != 3) return false;

            return center;
        }

        return false;
    }

    private void newIsomorphism(int first, int sec) {
        for (int i = 0; i < verticesNumber; i++){
            int temp = adj[first][i];
            adj[first][i] = adj[sec][i];
            adj[sec][i] = temp;
        }
        for (int i = 0; i < verticesNumber; i++){
            int temp = adj[i][first];
            adj[i][first] = adj[i][sec];
            adj[i][sec] = temp;
        }
    }

    public boolean isIsomorphic(Graph g){
        if (g.verticesNumber != verticesNumber || g.edgesNumber != edgesNumber) return false;
        if ((isSimple() && !g.isSimple()) || (g.isSimple() && !isSimple())) return false;
        if ((isConnected() && !g.isConnected()) || (g.isConnected() && !isConnected())) return false;
        if ((isCycle() && !g.isCycle()) || (g.isCycle() && !isCycle())) return false;

        if ((isCompleteBipartite() && !g.isCompleteBipartite()) ||
                (g.isCompleteBipartite() && !isCompleteBipartite())) return false;
        if(!compareDegrees(g)) return false;
        if (g.isWheel() && isWheel()) return true;
        if (g.isRegular() && isRegular()) return true;
        return checkMatrices(g, 0);
    }

    private boolean checkMatrices(Graph g, int k) {
        if (this.compareGraph(g)) return true;

        for (int i = k; i < verticesNumber - 1; i++)
            for (int j = i + 1; j < verticesNumber; j++){
                Graph temp = new Graph(Utils.copyMatrix(g.adj));
                temp.newIsomorphism(i, j);
                if (checkMatrices(temp, k + 1)) return true;
            }

        return false;
    }

    private boolean compareDegrees(Graph g) {
        HashMap<Integer, Integer> h = new HashMap<>();
        HashMap<Integer, Integer> h2 = new HashMap<>();

        for (int i = 0; i < verticesNumber; i++){
            int d = getVertexDegree(i);
            int d2 = g.getVertexDegree(i);
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

    public Set<Edge> edgeSet(){
        Set<Edge> edges = new HashSet<>();
        for (int i = 0; i < verticesNumber; i++)
            for (int j = 0; j < verticesNumber; j++)
                if (adj[i][j] > 0){
                    Edge e = new Edge(i, j);
                    edges.add(e);
                }
        return edges;
    }

    public boolean compareGraph(Graph g) {
        for (int i = 0; i < verticesNumber; i++)
            for (int j = 0; j < verticesNumber; j++)
                if (g.adj[i][j] != adj[i][j]) return false;
        return true;
    }

    private boolean isCycle() {
        if (verticesNumber == 1 && adj[0][0] == 2) return true;
        if (isDirected()){
            if (verticesNumber == 2 && edgesNumber == 2 && adj[0][1] == 1 && adj[1][0] == 1) return true;
            if (isSimple()){
                for (int i = 0; i < verticesNumber; i++)
                    if (getVertexDegree(i) != 1)
                        return false;
                return true;
            }
        }
        else{
            if (verticesNumber == 2 && edgesNumber == 2 && adj[0][1] == 2 && adj[1][0] == 2) return true;
            if (isSimple() && isConnected()) {
                for (int i = 0; i < verticesNumber; i++)
                    if (getVertexDegree(i) != 2)
                        return false;
                return true;
            }
        }
        return false;
    }

    private boolean isComplete() {
        if (isSimple()) {
            for (int i = 0; i < verticesNumber; i++)
                if (getVertexDegree(i) != verticesNumber - 1)
                    return false;
            return true;
        }
        return false;
    }

    private boolean isRegular() {
        if (isSimple()) {
            int deg = getVertexDegree(0);
            for (int i = 1; i < verticesNumber; i++)
                if (deg != getVertexDegree(i))
                    return false;
            return true;
        }
        return false;
    }

    public ArrayList<ArrayList<Integer>> getComponents(){
        ArrayList<ArrayList<Integer>> components = new ArrayList<>();
        boolean [] visited = new boolean[verticesNumber];
        int x = 0;
        boolean all = true;
        while (all) {
            ArrayList<Integer> component = new ArrayList<>();
            LinkedList<Integer> q = new LinkedList<>();
            visited[x] = true;
            q.add(x);
            while (q.size() != 0) {
                int c = q.poll();
                component.add(c);
                for (int i = 0; i < verticesNumber; i++) {
                    if (adj[c][i] > 0 && !visited[i]) {
                        visited[i] = true;
                        q.add(i);
                    }
                }
            }
            components.add(component);
            for (int i = 0; i < verticesNumber; i++)
                if (!visited[i]) {
                    all = true;
                    x = i;
                    break;
                }
                else
                    all = false;
        }
        return components;
    }

    private boolean containsCycle2(){
        if (containsLoop(adj)) return true;
        return getEdgesNumber() >= verticesNumber;
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
        return getRecursiveSpanningTree();
    }

    private int getRecursiveSpanningTree() {
        if (isConnected()){
            if (containsLoop()) {
                Graph temp = new Graph(Utils.copyMatrix(adj));
                temp.removeLoops();
                if (!temp.containsCycle2())
                    return 1;
            }
            if (!containsCycle2())
                return 1;
        }
        else return 0;

        int [] v = chooseVertices();
        Graph first = new Graph(Utils.copyMatrix(adj));
        Graph second = new Graph(Utils.copyMatrix(adj));
        first.edgeReduction(v[0], v[1]);
        second.edgeContraction(v[0], v[1]);

        return first.getRecursiveSpanningTree() + second.getRecursiveSpanningTree();
    }

    private void removeLoops() {
        for (int i = 0; i < adj.length; i++)
            adj[i][i] = 0;
    }

    private int[] chooseVertices() {
        Random r = new Random();
        int f, l;
        do {
            f = r.nextInt(verticesNumber);
            l = r.nextInt(verticesNumber);
        }while (adj[f][l] == 0 || f == l);
        return new int[]{f, l};
    }

    private void edgeReduction(int f, int l) {
        adj[f][l] = adj[f][l] - 1;
        adj[l][f] = adj[l][f] - 1;
    }

    private void edgeContraction(int f, int l){
        edgeReduction(f, l);
        for (int i = 0; i < verticesNumber; i++)
            adj[f][i] += adj[l][i];
        for (int i = 0; i < verticesNumber; i++)
            adj[i][f] += adj[i][l];

        int row = 0;
        for (int i = 0; i < verticesNumber; i++) {
            if (i != l) {
                int col = 0;
                for (int j = 0; j < verticesNumber; j++)
                    if (j != l) {
                        adj[row][col] = adj[i][j];
                        col++;
                    }
                row++;
            }
        }
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

    public int getVertexDegree(int vertex) {
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

    public ArrayList<Edge> eulerianWalk(){
        short numberOfOddVertices = 0;
        int chosenVertex = 0;
        for (int i = 0; i < verticesNumber; i++) {
            if (numberOfOddVertices > 2){
                System.out.println("not eulerian");
                return new ArrayList<>();
            }
            if (getVertexDegree(i) % 2 == 1) {
                chosenVertex = i;
                numberOfOddVertices++;
            }
        }

        Graph newGraph = new Graph(Utils.copyMatrix(adj));
        ArrayList<Edge> edgeSequence = new ArrayList<>();

        while (newGraph.edgesNumber > 0) {
            boolean cutEdge = true;
            int dest = -1;
            for (int i = 0; i < verticesNumber; i++) {
                Graph temp = new Graph(Utils.copyMatrix(newGraph.adj));
                if (temp.removeEdge(chosenVertex, i)){
                    if (temp.isConnected()){
                        newGraph = temp;
                        Edge edge = new Edge(chosenVertex, i);
                        edgeSequence.add(edge);
                        chosenVertex = i;
                        cutEdge = false;
                        break;
                    } else if (!temp.isConnected()) dest = i;
                }
            }
            if (cutEdge && dest != -1) {
                newGraph.removeEdge(chosenVertex, dest);
                Edge edge = new Edge(chosenVertex, dest);
                edgeSequence.add(edge);
                chosenVertex = dest;
            }
        }
        return edgeSequence;
    }

    public int maxMatchingNumber(){ //alpha prime
        if (isPath() || isCycle() || (verticesNumber % 2 == 0 && isComplete()))
            return verticesNumber / 2;
        if (isBipartite())
            if (isRegular())
                return verticesNumber / 2;
            else
                return minVertexCover();

        ArrayList<Integer> indices = new ArrayList<>();
        for (int i = 0; i < verticesNumber; i++)
            if (getVertexDegree(i) == 0) indices.add(i);
        if (!indices.isEmpty()) {
            Graph g = new Graph(Utils.copyMatrix(adj));
            g.removeSingleVertices(indices);
            return g.verticesNumber - g.minEdgeCover();
        }

        return verticesNumber - minEdgeCover();
    }

    private void removeSingleVertices(ArrayList<Integer> indices) {
        int l = verticesNumber - indices.size();
        int [][] a = new int[l][l];
        int count = 0;
        int counter = 0;
        for (int i = 0; i < verticesNumber; i++) {
            int count1 = 0;
            int counter1 = 0;
            if (i == indices.get(counter)){
                counter++;
                continue;
            }
            for (int j = 0; j < verticesNumber; j++)
                if (j == indices.get(counter1)) {
                    a[count][count1] = adj[i][j];
                    count1++;
                }
                else counter1++;

            count++;
        }
        adj = a;
        verticesNumber = adj.length;
    }

    public int minEdgeCover(){ //beta prime
        for (int i = 0; i < verticesNumber; i++)
            if (getVertexDegree(i) == 0) return 0;

        ArrayList<Edge> edges = new ArrayList<>(edgeSet());
        ArrayList<ArrayList<Edge>> subsets = Utils.subsets(edges, 2);
        for (ArrayList<Edge> es: subsets){
            boolean [] check = new boolean[verticesNumber];
            for (Edge e: es){
                check[e.getFirst()] = true;
                check[e.getDestination()] = true;
            }
            boolean checkAll = false;
            for (int i = 0; i < verticesNumber; i++) {
                if (!check[i]){
                    checkAll = false;
                    break;
                }
                checkAll = check[i];
            }
            if (checkAll) return es.size();
        }
        return 1;
    }

    public int maxIndependentSet(){ //alpha
        if (isComplete())
            return 1;
        if (isPath() || isCycle() || isCompleteBipartite())
            return verticesNumber / 2;

        ArrayList<ArrayList<Integer>> subsets = Utils.subsets(verticesNumber, 2);
        int max = 1;
        for (ArrayList<Integer> s : subsets){
            boolean isIndependent = true;
            for (int i = 0; i < s.size(); i++) {
                for (int j = i + 1; j < s.size(); j++)
                    if (adj[s.get(i)][s.get(j)] != 0) {
                        isIndependent = false;
                        break;
                    }
                if (!isIndependent) break;
            }
            if (isIndependent)
                if (s.size() > max) max = s.size();
        }

        return max;
    }

    public int minVertexCover(){ //beta
        return verticesNumber - maxIndependentSet();
    }

    @Override
    public String toString() {
        return "Graph{" +
                "adj=" + Arrays.toString(adj) +
                ", verticesNumber=" + verticesNumber +
                ", edgesNumber=" + edgesNumber +
                '}';
    }
}