import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.StringTokenizer;

public class Main {
    static BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
    static BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(System.out));
    static StringTokenizer st;
    private static List<double[]> position = new ArrayList<>();
    private static int height;
    private static int vertical;
    private static int xAxis;
    private static int yAxis;
    private static char[][] box;
    private static double[] white = new double[2];
    private static int cnt = 0;
    private static double result = 0;

    public static void main(String[] args) throws IOException {

        st = new StringTokenizer(br.readLine());
        height = Integer.parseInt(st.nextToken());
        vertical = Integer.parseInt(st.nextToken());
        xAxis = 0;
        yAxis = 0;
        box = new char[height][vertical];

        for (int i = 0; i < height; i++) {
            String str = br.readLine();
            for (int j = 0; j < vertical; j++) {
                box[i][j] = str.charAt(j);
                if (box[i][j] == '*') {
                    xAxis = j;
                    yAxis = i;
                    white = new double[]{j + 0.5, i + 0.5};
                } else if (box[i][j] == '#') {
                    cnt++;
                }
            }
        }

        List<double[]> containsed = new ArrayList<>();
        for (int i = 0; i <= 1; i++) {
            for (int j = 0; j <= 1; j++) {
                double x = vertical * j;
                double y = height * i;
                double[] item = new double[]{x - white[0] > 0 ? 1 : -1, (y - white[1]) / (x - white[0])};
                containsed.add(item);
            }
        }
        for (int i = 0; i <= height; i++) {
            for (int j = 0; j <= vertical; j++) {
                double[] item = new double[]{j - white[0] > 0 ? 1 : -1, (i - white[1]) / (j - white[0])};
                if (!containsed.contains(item)) {
                    containsed.add(item);
                }
            }
        }
        containsed.sort((a, b) -> a[0] != b[0] ? Double.compare(a[0], b[0]) : Double.compare(a[1], b[1]));
        containsed.forEach(i -> search(xAxis, yAxis, white, i));

        for (int i = 0; i < position.size(); i++) {
            int j = (i + 1) % (position.size());
            double a = white[0] * position.get(i)[1] + position.get(i)[0] * position.get(j)[1] + position.get(j)[0] * white[1];
            double b = position.get(i)[0] * white[1] + position.get(j)[0] * position.get(i)[1] + white[0] * position.get(j)[1];
            result += Math.abs(a - b) / 2;
        }
        System.out.printf("%.13f", vertical * height - result - cnt);
    }

    private static boolean height(int x, int y) {
        return y == -1 || y == height || box[y][x] == '#';
    }

    private static boolean vertical(int x, int y) {
        return x == -1 || x == vertical || box[y][x] == '#';
    }

    private static void search(int x, int y, double[] pos, double[] a) {
        if (x == -1 || x == vertical || y == -1 || y == height || box[y][x] == '#') {
            position.add(pos);
        } else {
            int x1 = a[0] > 0 ? x + 1 : x - 1;
            int x2 = a[0] > 0 ? x + 1 : x;
            int y1 = a[0] * a[1] > 0 ? y + 1 : y - 1;
            int y2 = a[0] * a[1] > 0 ? y + 1 : y;
            double y_nx = (pos[1] + a[1] * (x2 - pos[0]));
            double x_ny = (pos[0] + ((y2 - pos[1]) / a[1]));
            if (y + 10e-7 < y_nx && y_nx < y + 1 - 10e-7 && (x_ny < x - 10e-7 || x_ny > x + 1 + 10e-7)) {
                pos = new double[]{x2, y_nx};
                search(x1, y, pos, a);
            } else if (x + 10e-7 < x_ny && x_ny < x + 1 - 10e-7 && (y - 10e-7 > y_nx || y_nx > y + 1 + 10e-7)) {
                pos = new double[]{x_ny, y2};
                search(x, y1, pos, a);
            } else {
                pos = new double[]{x2, y2};
                if (a[1] > 0 && vertical(x1, y)) position.add(pos);
                if (a[1] < 0 && height(x, y1)) position.add(pos);
                search(x1, y1, pos, a);
                if (a[1] > 0 && height(x, y1)) position.add(pos);
                if (a[1] < 0 && vertical(x1, y)) position.add(pos);
            }
        }
    }
}
