import java.util.Objects;

public class Edge {
    private final int first;
    private final int destination;

    public Edge(int first, int destination){
        this.first = first;
        this.destination = destination;
    }

    public int getDestination() {
        return destination;
    }

    public int getFirst() {
        return first;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Edge edge = (Edge) o;
        return (first == edge.first && destination == edge.destination) ||
                (first == edge.destination && destination == edge.first);
    }

    @Override
    public int hashCode() {
        return 1;
    }

    @Override
    public String toString() {
        return "{" +
                "" + first +
                "," + destination +
                '}';
    }
}
