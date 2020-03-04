import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class Solver {

    public static final String ANSI_BLACK_BACKGROUND = "\u001B[40m";
    public static final String ANSI_RED_BACKGROUND = "\u001B[41m";
    public static final String ANSI_GREEN_BACKGROUND = "\u001B[42m";
    public static final String ANSI_YELLOW_BACKGROUND = "\u001B[43m";
    public static final String ANSI_BLUE_BACKGROUND = "\u001B[44m";
    public static final String ANSI_PURPLE_BACKGROUND = "\u001B[45m";
    public static final String ANSI_CYAN_BACKGROUND = "\u001B[46m";
    public static final String ANSI_WHITE_BACKGROUND = "\u001B[47m";
    public static final String ANSI_POS_BACKGROUND = "\033[48;5;23m";


    public static final String ANSI_RESET = "\u001B[0m";
    public static final String ANSI_BLACK = "\u001B[30m";
    public static final String ANSI_RED = "\u001B[31m";
    public static final String ANSI_RED2 = "\033[38;5;213m";
    public static final String ANSI_GREEN = "\033[38;5;28m";
    public static final String ANSI_YELLOW = "\u001B[33m";
    public static final String ANSI_BLUE = "\u001B[34m";
    public static final String ANSI_PURPLE = "\u001B[35m";
    public static final String ANSI_CYAN = "\u001B[36m";
    public static final String ANSI_WHITE = "\u001B[37m";



    /** ------------ Параметры среды --------------*/
    private static final int r = 100; // радиус
    private static final int h = 100; // высота
    private static final double rho = 7900.0; // плотность
    private static final double lamda = 1.0; // параметр Ляме
    private static final double mu = 7.89e+10; // модуль сдвига (Тоже параметр Ляме)
    private static final double v_medium = 1.0; // скорость звука в среде
    private static final double y_0 = 245e+4; // предел текучести

    /** ------------ Параметры удара --------------*/
    private static final int r_0 = 4;
    private static final double v_0 = -500.0;

    /** ------------ Параметры солверов --------------*/
    private static final double viscosity_coff = 0.5;
    private static final double timestep = 0.6666666666e-7;
    private static final double CFL = 0.5;
    private static final int x_count = 11;
    private static final int y_count = 22;
    private static final double l_x = 0.1;
    private static final double l_y = 0.1;

//    private static double node[] = new double[x_res];
//    private static double cell[] = new double[x_res];
//    private static ArrayList<HashMap> node = Arrays.asList(new HashMap<String, Double>[h_res]);
//    private static List<Map<String , Double>> nodes  = new ArrayList<>();
//    private static List<Double> cells = Arrays.asList(new Double[h_res]);
    private static final Node[][] nodes = new Node[x_count][y_count];
    private static Node[][] oldNodes = new Node[x_count][y_count];
    private static Cell[][] cells = new Cell[x_count][y_count];
    private static Cell[][] oldCells = new Cell[x_count][y_count];

    public Solver() throws WrongAccesException, WrongSpeedException, IOException {
        System.out.println(ANSI_RED+"-------Solver started--------"+ANSI_RESET);

        System.out.println(ANSI_RED2+"======== Prepare step started ========"+ANSI_RESET);
        startPrepStep();

        System.out.println(ANSI_RED2+"======== Solver step started ========"+ANSI_RESET);
        startSolveStep();
    }

    public static void main(String[] args) throws WrongAccesException, WrongSpeedException, IOException {
        new Solver();
    }

    private void writeToFile() throws IOException {
        File myObj = new File("./pts.txt");
        File myObj1 = new File("./prims.txt");
        FileWriter myWriter = new FileWriter("./pts.txt");
        StringBuilder f = new StringBuilder();
        f.append(x_count).append("_").append(y_count).append("\n");
        for(int k=0; k<x_count; k++)//проходим все слои, выставляя координаты(прямоугольная недеформированная сетка), скорость ставим ноль, записываем номер
        {
            for (int j = 0; j < y_count; j++) {
                f.append("pt_").append(Arrays.toString(nodes[k][j].getNum())).append("_").append(Arrays.toString(nodes[k][j].getvCoords())).append("_").append(Arrays.toString(nodes[k][j].getvV())).append("\n");
            }
        }
        myWriter.write(f.toString());
        myWriter.close();

        myWriter = new FileWriter("./prims.txt");
        f = new StringBuilder();
//        f.append(x_count).append("_").append(y_count).append("\n");
        for(int k=0; k<x_count; k++)//проходим все слои, выставляя координаты(прямоугольная недеформированная сетка), скорость ставим ноль, записываем номер
        {
            for (int j = 0; j < y_count; j++) {
                f.append("prim_").append(Arrays.toString(cells[k][j].getStress())).append("\n");
            }
        }
        myWriter.write(f.toString());
        myWriter.close();
    }

    private void startPrepStep() throws WrongAccesException {
        /** Записываем начальные условия для узлов
         * В узле хранятся:
         * - координтаы
         * - скорость
         * - phi, влияет на скорости
         * - alpha, влияет на скорости
         * - beta, влияет на скорости
         * */

        double x_cell_length = l_x/(x_count-2);
        double y_cell_length = l_y/(y_count-2);

        /**
         * Заполняем нулями
         * */
        double m=-1; //итератор по Х
        double i=-1; // итератор по Y

        System.out.println(ANSI_RED+"---->> Filling grid with zeroes"+ANSI_RESET);

        for(int k=0; k<x_count; k++)//проходим все слои, выставляя координаты(прямоугольная недеформированная сетка), скорость ставим ноль, записываем номер
        {
            for(int j=0; j<y_count; j++)
            {
                nodes[k][j] = new Node(new double[]{m*x_cell_length, i*y_cell_length}, new double[]{v_0, 0}, new int[]{k,j});
                oldNodes[k][j] = new Node(new double[]{m*x_cell_length, i*y_cell_length}, new double[]{v_0, 0}, new int[]{k,j});

                /**
                 *  flags: 0 - принадлежит телу, 1 - мнимая, 2 - отраженная, 3 4 5 6 - углы, 7 8 9 10 - стенки
                 *
                 *  1 4 9 9 9 9 9 9 9 9 9 5
                 *  1 8 0 0 0 0 0 0 0 0 0 10
                 *  1 8 0 0 0 0 0 0 0 0 0 10
                 *  1 8 0 0 0 0 0 0 0 0 0 10
                 *  1 8 0 0 0 0 0 0 0 0 0 10
                 *  1 3 7 7 7 7 7 7 7 7 7 6
                 *  1 2 2 2 2 2 2 2 2 2 2 2
                 */

                if(k == 0){
                    nodes[k][j].setFlag(1);
                    oldNodes[k][j].setFlag(1);
                }else if(j == 0){
                    nodes[k][j].setFlag(2);
                    oldNodes[k][j].setFlag(2);
                }else if(k == 1 && j == 1){
                    nodes[k][j].setFlag(3);
                    oldNodes[k][j].setFlag(3);
                }else if(k == 1 && j == y_count-1){
                    nodes[k][j].setFlag(4);
                    oldNodes[k][j].setFlag(4);
                }else if(k == x_count-1 && j == y_count-1){
                    nodes[k][j].setFlag(5);
                    oldNodes[k][j].setFlag(5);
                }else if(k == x_count-1  && j == 1){
                    nodes[k][j].setFlag(6);
                    oldNodes[k][j].setFlag(6);
                }else if(k == 1){
                    nodes[k][j].setFlag(8);
                    oldNodes[k][j].setFlag(8);
                }else if(j == 1){
                    nodes[k][j].setFlag(7);
                    oldNodes[k][j].setFlag(7);
                }else if(j == y_count-1){
                    nodes[k][j].setFlag(9);
                    oldNodes[k][j].setFlag(9);
                }else if(k == x_count-1){
                    nodes[k][j].setFlag(10);
                    oldNodes[k][j].setFlag(10);
                }
                i++;
            }
            m++;
            i=-1;
        }

        System.out.println(ANSI_GREEN+"---->> OK"+ANSI_RESET);

        printFields(false, false, true, false, false);

        /**
         *Задаем скорости в момент соприкосновения
         */
        System.out.println(ANSI_RED+"---->> Applying initial hit values"+ANSI_RESET);

        for(int iter = 0; iter < y_count; iter++){
            nodes[0][iter].setV(0.0, 0.0);
            oldNodes[0][iter].setV(0.0, 0.0);
            nodes[1][iter].setV(0.0, 0.0);
            oldNodes[1][iter].setV(0.0, 0.0);
//            nodes[iter][0].setV(v_0, 0.0); /** отражение по Y*/
        }
        System.out.println(ANSI_GREEN+"---->> OK"+ANSI_RESET);

        /**
         *Создаем ячейки
         */
        double mass = (2.0/3.0)*rho*Math.PI;
        for(int k=0; k<x_count; k++)//проходим все слои, выставляя координаты(прямоугольная недеформированная сетка), скорость ставим ноль, записываем номер
        {
            for(int j=0; j<y_count; j++)
            {
                cells[k][j] = new Cell(new int[]{k,j});
                oldCells[k][j] = new Cell(new int[]{k,j});

                if(k!=0 && j!=0){
                    setAreaAndMass(k, j, mass);
                }

                if(k == x_count-1 || j == y_count-1){
                    cells[k][j].setFlag(3);
                    oldCells[k][j].setFlag(3);
                } else if(k == 0) {
                    cells[k][j].setFlag(1);
                    oldCells[k][j].setFlag(1);
                } else if(j == 0){
                    cells[k][j].setFlag(2);
                    oldCells[k][j].setFlag(2);
                }
            }
        }

        /**
         *У отраженных ячеек перезаписываем массу
         */
        for(int it = 1; it < x_count-1; it++){
            cells[it][0].setMass(cells[it][1].getMass());
        }
    }

    private void startSolveStep() throws WrongAccesException, WrongSpeedException, IOException {
        System.out.println(ANSI_BLUE+"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"+ANSI_RESET);
        System.out.println(ANSI_BLUE+"~~~~~~~~~~~~~~~~~~ Initial valuse: ~~~~~~~~~~~~~~~~~~"+ANSI_RESET);
        System.out.println(ANSI_BLUE+"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"+ANSI_RESET);
        printFields(true, true, false, false, false);

        copyTimeStep();

        double phi;
        double alpha;
        double betta;

        System.out.println(ANSI_BLUE+"\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"+ANSI_RESET);
        System.out.println(ANSI_BLUE+"~~~~~~~~~~~~~~~~~~ 1-st step valuse: ~~~~~~~~~~~~~~~~~~"+ANSI_RESET);
        System.out.println(ANSI_BLUE+"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"+ANSI_RESET);

        printFields(false, false, false, false, true);

        for(int k=1; k < x_count; k++) {
            for (int j = 1; j < y_count; j++) {
                if (oldNodes[k][j].getFlag() == 9) {
                    /**
                     * TODO: Разобраться со множителями phi, alpha, betta
                     */
                    phi = 0.25 * (
                            oldCells[k - 1][j - 1].getArea() * oldCells[k - 1][j - 1].getDensity() / oldCells[k - 1][j - 1].getRelDensity() +
                                    oldCells[k][j - 1].getArea() * oldCells[k][j - 1].getDensity() / oldCells[k][j - 1].getRelDensity()
                    );

                    alpha = 0.5 * (
                            oldCells[k - 1][j - 1].getArea() * oldCells[k - 1][j - 1].getStress()[2] / oldCells[k - 1][j - 1].getMass() +
                                    oldCells[k][j - 1].getArea() * oldCells[k][j - 1].getStress()[2] / oldCells[k][j - 1].getMass()
                    );

                    betta = 0.5 * (
                            oldCells[k - 1][j - 1].getArea() * (oldCells[k - 1][j - 1].getStress()[1] - oldCells[k - 1][j - 1].getStress()[3]) / oldCells[k - 1][j - 1].getMass() +
                                    oldCells[k][j - 1].getArea() * (oldCells[k][j - 1].getStress()[1] - oldCells[k][j - 1].getStress()[3]) / oldCells[k][j - 1].getMass()
                    );

                    double xv = oldNodes[k][j].getV()[0] - (timestep / (2 * phi)) * (
                            oldCells[k - 1][j - 1].getStress()[0] * (oldNodes[k - 1][j].getCoords()[1] - oldNodes[k][j - 1].getCoords()[1]) +
                                    oldCells[k][j - 1].getStress()[0] * (oldNodes[k][j - 1].getCoords()[1] - oldNodes[k + 1][j].getCoords()[1]) -

                                    oldCells[k - 1][j - 1].getStress()[2] * (oldNodes[k - 1][j].getCoords()[0] - oldNodes[k][j - 1].getCoords()[0]) -
                                    oldCells[k][j - 1].getStress()[2] * (oldNodes[k][j - 1].getCoords()[0] - oldNodes[k + 1][j].getCoords()[0])) +
                            timestep * alpha;

                    double yv = oldNodes[k][j].getV()[1] + (timestep / (2 * phi)) * (
                            oldCells[k - 1][j - 1].getStress()[1] * (oldNodes[k - 1][j].getCoords()[0] - oldNodes[k][j - 1].getCoords()[0]) +
                                    oldCells[k][j - 1].getStress()[1] * (oldNodes[k][j - 1].getCoords()[0] - oldNodes[k + 1][j].getCoords()[0]) -

                                    oldCells[k - 1][j - 1].getStress()[2] * (oldNodes[k - 1][j].getCoords()[0] - oldNodes[k][j - 1].getCoords()[0]) -
                                    oldCells[k][j - 1].getStress()[2] * (oldNodes[k][j - 1].getCoords()[0] - oldNodes[k + 1][j].getCoords()[0])) +
                            timestep * betta;
                    nodes[k][j].setV(xv, yv);
                    nodes[k][j].setCoords(oldNodes[k][j].getCoords()[0] + xv * timestep, oldNodes[k][j].getCoords()[1] + yv * timestep);
                    nodes[k][j].setColor(nodes[k][j].getColor() + "\u001b[4m");

                } else if (oldNodes[k][j].getFlag() == 10) {
                    /**
                     * TODO: Разобраться со множителями phi, alpha, betta
                     */
                    phi = 0.25 * (
                            oldCells[k - 1][j].getArea() * oldCells[k - 1][j].getDensity() / oldCells[k - 1][j].getRelDensity() +
                                    oldCells[k - 1][j - 1].getArea() * oldCells[k - 1][j - 1].getDensity() / oldCells[k - 1][j - 1].getRelDensity()
                    );

                    alpha = 0.5 * (
                            oldCells[k - 1][j].getArea() * oldCells[k - 1][j].getStress()[2] / oldCells[k - 1][j].getMass() +
                                    oldCells[k - 1][j - 1].getArea() * oldCells[k - 1][j - 1].getStress()[2] / oldCells[k - 1][j - 1].getMass()
                    );

                    if (Double.isNaN(alpha)) {
                        System.out.println("ERROR!!!!!!!1");
                        System.out.println(oldCells[k][j - 1].getArea());
                        System.out.println(oldCells[k][j - 1].getStress()[2]);
                        System.out.println(oldCells[k][j - 1].getMass());
                        System.out.println(k + "  " + (j - 1));
                        throw new WrongSpeedException();
                    }

                    betta = 0.5 * (
                            oldCells[k - 1][j].getArea() * (oldCells[k - 1][j].getStress()[1] - oldCells[k - 1][j].getStress()[3]) / oldCells[k - 1][j].getMass() +
                                    oldCells[k - 1][j - 1].getArea() * (oldCells[k - 1][j - 1].getStress()[1] - oldCells[k - 1][j - 1].getStress()[3]) / oldCells[k - 1][j - 1].getMass()
                    );

                    double xv = oldNodes[k][j].getV()[0] - (timestep / (2 * phi)) * (
                            oldCells[k - 1][j].getStress()[0] * (oldNodes[k][j + 1].getCoords()[1] - oldNodes[k - 1][j].getCoords()[1]) +
                                    oldCells[k - 1][j - 1].getStress()[0] * (oldNodes[k - 1][j].getCoords()[1] - oldNodes[k][j - 1].getCoords()[1]) -

                                    oldCells[k - 1][j].getStress()[2] * (oldNodes[k][j + 1].getCoords()[0] - oldNodes[k - 1][j].getCoords()[0]) -
                                    oldCells[k - 1][j - 1].getStress()[2] * (oldNodes[k - 1][j].getCoords()[0] - oldNodes[k][j - 1].getCoords()[0])
                    ) +
                            timestep * alpha;

                    double yv = oldNodes[k][j].getV()[1] + (timestep / (2 * phi)) * (
                            oldCells[k - 1][j].getStress()[1] * (oldNodes[k][j + 1].getCoords()[0] - oldNodes[k - 1][j].getCoords()[0]) +
                                    oldCells[k - 1][j - 1].getStress()[1] * (oldNodes[k - 1][j].getCoords()[0] - oldNodes[k][j - 1].getCoords()[0]) +

                                    oldCells[k - 1][j].getStress()[2] * (oldNodes[k][j + 1].getCoords()[0] - oldNodes[k - 1][j].getCoords()[0]) -
                                    oldCells[k - 1][j - 1].getStress()[2] * (oldNodes[k - 1][j].getCoords()[0] - oldNodes[k][j - 1].getCoords()[0])
                    ) +
                            timestep * betta;
                    nodes[k][j].setV(xv, yv);
                    nodes[k][j].setCoords(oldNodes[k][j].getCoords()[0] + xv * timestep, oldNodes[k][j].getCoords()[1] + yv * timestep);
                    nodes[k][j].setColor(nodes[k][j].getColor() + "\u001b[4m");

                } else if (nodes[k][j].getFlag() == 4) {
                    phi = 0.25 * (
                            oldCells[k][j - 1].getArea() * oldCells[k][j - 1].getDensity() / oldCells[k][j - 1].getRelDensity()
                    );


                    betta = (
                            oldCells[k][j - 1].getArea() * (oldCells[k][j - 1].getStress()[1] - oldCells[k][j - 1].getStress()[3]) / oldCells[k][j - 1].getMass()
                    );

                    if (Double.isNaN(betta)) {
                        System.out.println("ERROR!!!!!!!1");
                        System.out.println(oldCells[k][j - 1].getArea());
                        System.out.println(oldCells[k][j - 1].getStress()[3]);
                        System.out.println(oldCells[k - 1][j].getMass());
                        System.out.println(k + "  " + (j - 1));
                        throw new WrongSpeedException(String.valueOf(nodes[k][j].getFlag()));
                    }

                    double yv = oldNodes[k][j].getV()[1] + (timestep / (2 * phi)) * (
                            oldCells[k - 1][j - 1].getStress()[1] * (oldNodes[k - 1][j].getCoords()[0] - oldNodes[k][j - 1].getCoords()[0]) +
                                    oldCells[k][j - 1].getStress()[1] * (oldNodes[k][j - 1].getCoords()[0] - oldNodes[k + 1][j].getCoords()[0]) -

                                    oldCells[k - 1][j - 1].getStress()[2] * (oldNodes[k - 1][j].getCoords()[0] - oldNodes[k][j - 1].getCoords()[0]) -
                                    oldCells[k][j - 1].getStress()[2] * (oldNodes[k][j - 1].getCoords()[0] - oldNodes[k + 1][j].getCoords()[0])) +
                            timestep * betta;
                    nodes[k][j].setV(0, yv);
                    nodes[k][j].setCoords(oldNodes[k][j].getCoords()[0], oldNodes[k][j].getCoords()[1] + yv * timestep);
                    nodes[k][j].setColor(nodes[k][j].getColor() + "\u001b[4m");

                } else if (nodes[k][j].getFlag() == 6) {
                    phi = 0.25 * (
                            oldCells[k - 1][j].getArea() * oldCells[k - 1][j].getDensity() / oldCells[k - 1][j].getRelDensity()
                    );

                    alpha = (
                            oldCells[k - 1][j].getArea() * oldCells[k - 1][j].getStress()[2] / oldCells[k - 1][j].getMass()
                    );

                    if (Double.isNaN(alpha)) {
                        System.out.println("ERROR!!!!!!!1");
                        System.out.println(oldCells[k][j - 1].getArea());
                        System.out.println(oldCells[k][j - 1].getStress()[2]);
                        System.out.println(oldCells[k][j - 1].getMass());
                        System.out.println(k + "  " + (j - 1));
                        throw new WrongSpeedException(String.valueOf(oldCells[k][j].getFlag()));
                    }

                    if (alpha != 0) {
                        throw new WrongSpeedException(String.valueOf(oldCells[k][j].getFlag()));
                    }

                    double xv = oldNodes[k][j].getV()[0] - (timestep / (2 * phi)) * (
                            oldCells[k - 1][j].getStress()[0] * (oldNodes[k][j + 1].getCoords()[1] - oldNodes[k - 1][j].getCoords()[1]) +
                                    oldCells[k - 1][j - 1].getStress()[0] * (oldNodes[k - 1][j].getCoords()[1] - oldNodes[k][j - 1].getCoords()[1]) -

                                    oldCells[k - 1][j].getStress()[2] * (oldNodes[k][j + 1].getCoords()[0] - oldNodes[k - 1][j].getCoords()[0]) -
                                    oldCells[k - 1][j - 1].getStress()[2] * (oldNodes[k - 1][j].getCoords()[0] - oldNodes[k][j - 1].getCoords()[0])
                    ) +
                            timestep * alpha;

                    if (Double.isNaN(xv)) {
                        throw new WrongSpeedException();
                    }

                    nodes[k][j].setV(xv, 0);
                    nodes[k][j].setCoords(oldNodes[k][j].getCoords()[0] + xv * timestep, oldNodes[k][j].getCoords()[1]);
                    nodes[k][j].setColor(nodes[k][j].getColor() + "\u001b[4m");

                } else if (nodes[k][j].getFlag() == 5) {
                    phi = 0.25 * (
                            oldCells[k - 1][j - 1].getArea() * oldCells[k - 1][j - 1].getDensity() / oldCells[k - 1][j - 1].getRelDensity()
                    );

                    alpha = 0.25 * (
                            oldCells[k - 1][j - 1].getArea() * oldCells[k - 1][j - 1].getStress()[2] / oldCells[k - 1][j - 1].getMass()
                    );

                    if (Double.isNaN(alpha)) {
                        System.out.println("ERROR!!!!!!!1");
                        System.out.println(oldCells[k][j - 1].getArea());
                        System.out.println(oldCells[k][j - 1].getStress()[2]);
                        System.out.println(oldCells[k][j - 1].getMass());
                        System.out.println(k + "  " + (j - 1));
                        throw new WrongSpeedException(String.valueOf(oldCells[k][j].getFlag()));
                    }

                    betta = 0.25 * (
                            oldCells[k - 1][j - 1].getArea() * (oldCells[k - 1][j - 1].getStress()[1] - oldCells[k - 1][j - 1].getStress()[3]) / oldCells[k - 1][j - 1].getMass()
                    );

                    double xv = oldNodes[k][j].getV()[0] - (timestep / (2 * phi)) * (
                            oldCells[k - 1][j - 1].getStress()[0] * (oldNodes[k - 1][j].getCoords()[1] - oldNodes[k][j - 1].getCoords()[1]) -

                                    oldCells[k - 1][j - 1].getStress()[2] * (oldNodes[k - 1][j].getCoords()[0] - oldNodes[k][j - 1].getCoords()[0])) +
                            timestep * alpha;

                    double yv = oldNodes[k][j].getV()[1] + (timestep / (2 * phi)) * (
                            oldCells[k - 1][j - 1].getStress()[1] * (oldNodes[k - 1][j].getCoords()[0] - oldNodes[k][j - 1].getCoords()[0]) -

                                    oldCells[k - 1][j - 1].getStress()[2] * (oldNodes[k - 1][j].getCoords()[0] - oldNodes[k][j - 1].getCoords()[0])
                    ) +
                            timestep * betta;
                    nodes[k][j].setV(xv, yv);
                    nodes[k][j].setCoords(oldNodes[k][j].getCoords()[0] + xv * timestep, oldNodes[k][j].getCoords()[1] + yv * timestep);
                    nodes[k][j].setColor(nodes[k][j].getColor() + "\u001b[4m");

                } else {
                    phi = 0.25 * (
                            oldCells[k][j].getArea() * oldCells[k][j].getDensity() / oldCells[k][j].getRelDensity() +
                                    oldCells[k - 1][j].getArea() * oldCells[k - 1][j].getDensity() / oldCells[k - 1][j].getRelDensity() +
                                    oldCells[k - 1][j - 1].getArea() * oldCells[k - 1][j - 1].getDensity() / oldCells[k - 1][j - 1].getRelDensity() +
                                    oldCells[k][j - 1].getArea() * oldCells[k][j - 1].getDensity() / oldCells[k][j - 1].getRelDensity()
                    );

                    alpha = 0.25 * (
                            oldCells[k][j].getArea() * oldCells[k][j].getStress()[2] / oldCells[k][j].getMass() +
                                    oldCells[k - 1][j].getArea() * oldCells[k - 1][j].getStress()[2] / oldCells[k - 1][j].getMass() +
                                    oldCells[k - 1][j - 1].getArea() * oldCells[k - 1][j - 1].getStress()[2] / oldCells[k - 1][j - 1].getMass() +
                                    oldCells[k][j - 1].getArea() * oldCells[k][j - 1].getStress()[2] / oldCells[k][j - 1].getMass()
                    );

                    if (Double.isNaN(alpha)) {
                        System.out.println("ERROR!!!!!!!1");
                        System.out.println(oldCells[k][j - 1].getArea());
                        System.out.println(oldCells[k][j - 1].getStress()[2]);
                        System.out.println(oldCells[k][j - 1].getMass());
                        System.out.println(k + "  " + (j - 1));
                        throw new WrongSpeedException(String.valueOf(oldCells[k][j].getFlag()));
                    }

                    betta = 0.25 * (
                            oldCells[k][j].getArea() * (oldCells[k][j].getStress()[1] - oldCells[k][j].getStress()[3]) / oldCells[k][j].getMass() +
                                    oldCells[k - 1][j].getArea() * (oldCells[k - 1][j].getStress()[1] - oldCells[k - 1][j].getStress()[3]) / oldCells[k - 1][j].getMass() +
                                    oldCells[k - 1][j - 1].getArea() * (oldCells[k - 1][j - 1].getStress()[1] - oldCells[k - 1][j - 1].getStress()[3]) / oldCells[k - 1][j - 1].getMass() +
                                    oldCells[k][j - 1].getArea() * (oldCells[k][j - 1].getStress()[1] - oldCells[k][j - 1].getStress()[3]) / oldCells[k][j - 1].getMass()
                    );

                    double xv = oldNodes[k][j].getV()[0] - (timestep / (2 * phi)) * (oldCells[k][j].getStress()[0] * (oldNodes[k + 1][j].getCoords()[1] - oldNodes[k][j + 1].getCoords()[1]) +
                            oldCells[k - 1][j].getStress()[0] * (oldNodes[k][j + 1].getCoords()[1] - oldNodes[k - 1][j].getCoords()[1]) +
                            oldCells[k - 1][j - 1].getStress()[0] * (oldNodes[k - 1][j].getCoords()[1] - oldNodes[k][j - 1].getCoords()[1]) +
                            oldCells[k][j - 1].getStress()[0] * (oldNodes[k][j - 1].getCoords()[1] - oldNodes[k + 1][j].getCoords()[1]) -

                            oldCells[k][j].getStress()[2] * (oldNodes[k + 1][j].getCoords()[0] - oldNodes[k][j + 1].getCoords()[0]) -
                            oldCells[k - 1][j].getStress()[2] * (oldNodes[k][j + 1].getCoords()[0] - oldNodes[k - 1][j].getCoords()[0]) -
                            oldCells[k - 1][j - 1].getStress()[2] * (oldNodes[k - 1][j].getCoords()[0] - oldNodes[k][j - 1].getCoords()[0]) -
                            oldCells[k][j - 1].getStress()[2] * (oldNodes[k][j - 1].getCoords()[0] - oldNodes[k + 1][j].getCoords()[0])) +
                            timestep * alpha;

                    double yv = oldNodes[k][j].getV()[1] + (timestep / (2 * phi)) * (oldCells[k][j].getStress()[1] * (oldNodes[k + 1][j].getCoords()[0] - oldNodes[k][j + 1].getCoords()[0]) +
                            oldCells[k - 1][j].getStress()[1] * (oldNodes[k][j + 1].getCoords()[0] - oldNodes[k - 1][j].getCoords()[0]) +
                            oldCells[k - 1][j - 1].getStress()[1] * (oldNodes[k - 1][j].getCoords()[0] - oldNodes[k][j - 1].getCoords()[0]) +
                            oldCells[k][j - 1].getStress()[1] * (oldNodes[k][j - 1].getCoords()[0] - oldNodes[k + 1][j].getCoords()[0]) -

                            oldCells[k][j].getStress()[2] * (oldNodes[k + 1][j].getCoords()[0] - oldNodes[k][j + 1].getCoords()[0]) -
                            oldCells[k - 1][j].getStress()[2] * (oldNodes[k][j + 1].getCoords()[0] - oldNodes[k - 1][j].getCoords()[0]) -
                            oldCells[k - 1][j - 1].getStress()[2] * (oldNodes[k - 1][j].getCoords()[0] - oldNodes[k][j - 1].getCoords()[0]) -
                            oldCells[k][j - 1].getStress()[2] * (oldNodes[k][j - 1].getCoords()[0] - oldNodes[k + 1][j].getCoords()[0])) +
                            timestep * betta;
                    nodes[k][j].setV(xv, yv);
                    nodes[k][j].setCoords(oldNodes[k][j].getCoords()[0] + xv * timestep, oldNodes[k][j].getCoords()[1] + yv * timestep);
                    nodes[k][j].setColor(nodes[k][j].getColor() + "\u001b[4m");

                }

                if (k < x_count - 1 && j < y_count - 1) {
                    /**
                     * Deformation part
                     */

                    double A_a = (nodes[k + 1][j].getCoords()[0] * (nodes[k + 1][j + 1].getCoords()[1] - nodes[k][j + 1].getCoords()[1]) +
                            nodes[k + 1][j + 1].getCoords()[0] * (nodes[k][j + 1].getCoords()[1] - nodes[k + 1][j].getCoords()[1]) +
                            nodes[k][j + 1].getCoords()[0] * (nodes[k + 1][j].getCoords()[1] - nodes[k + 1][j + 1].getCoords()[1])) * 0.5;

                    double A_b = (nodes[k + 1][j].getCoords()[0] * (nodes[k][j + 1].getCoords()[1] - nodes[k][j].getCoords()[1]) +
                            nodes[k][j + 1].getCoords()[0] * (nodes[k][j].getCoords()[1] - nodes[k + 1][j].getCoords()[1]) +
                            nodes[k][j].getCoords()[0] * (nodes[k + 1][j].getCoords()[1] - nodes[k][j + 1].getCoords()[1])) * 0.5;

                    double V_tmp = (nodes[k + 1][j].getCoords()[1] + nodes[k + 1][j + 1].getCoords()[1] + nodes[k][j + 1].getCoords()[1]) *
                            A_a + (nodes[k + 1][j].getCoords()[1] + nodes[k][j].getCoords()[1] + nodes[k][j + 1].getCoords()[1]) *
                            A_b;

                    cells[k][j].setArea(A_a + A_b);
                    cells[k][j].setDensity((3.0 * oldCells[k][j].getMass()) / (2.0 * V_tmp));

                    double hlafstepA = (oldCells[k][j].getArea() + A_a + A_b) * 0.5;
                    double defXX = op1(0, k, j) / (2.0 * hlafstepA);
                    double defYY = op1(3, k, j) / (2.0 * hlafstepA);

                    double defTT = (cells[k][j].getRelDensity() - oldCells[k][j].getRelDensity()) / ((cells[k][j].getRelDensity() + oldCells[k][j].getRelDensity()) * 0.5 * timestep) - (defXX + defYY);

                    double defXY = (op1(1, k, j) + op1(2, k, j)) / (2.0 * hlafstepA);

                    defTT *= timestep;
                    defXY *= timestep;
                    defXX *= timestep;
                    defYY *= timestep;

                    double tmpVar = 2 * ((cells[k][j].getRelDensity() - oldCells[k][j].getRelDensity()) / (3.0 * (cells[k][j].getRelDensity() + oldCells[k][j].getRelDensity())));

                    cells[k][j].setDev(0,
                            oldCells[k][j].getDev()[0] + 2.0 * mu * (defXX - tmpVar)
                    );

                    cells[k][j].setDev(1,
                            oldCells[k][j].getDev()[1] + 2.0 * mu * (defYY - tmpVar)
                    );

                    cells[k][j].setDev(3,
                            oldCells[k][j].getDev()[3] + 2.0 * mu * (defTT - tmpVar)
                    );

                    cells[k][j].setStress(2, oldCells[k][j].getStress()[2] + mu * defXY);

                    double fnext = 2*Math.pow(cells[k][j].getStress()[2], 2) + Math.pow(cells[k][j].getDev()[1], 2) + Math.pow(cells[k][j].getDev()[3], 2) + Math.pow(cells[k][j].getDev()[0], 2);

                    if(fnext > Math.pow(y_0, 2)*2.0/3.0){
                        double mul = y_0*(Math.sqrt(2.0/3.0)*fnext);
                        cells[k][j].setStress(2, cells[k][j].getStress()[2]*mul);
                        cells[k][j].setDev(0, cells[k][j].getDev()[0]*mul);
                        cells[k][j].setDev(1, cells[k][j].getDev()[1]*mul);
                        cells[k][j].setDev(3, cells[k][j].getDev()[3]*mul);
                    }
                    cells[k][j].setStress(0,cells[k][j].getDev()[0] + cells[k][j].getPressure());
                    cells[k][j].setStress(1,cells[k][j].getDev()[1] + cells[k][j].getPressure());
                    cells[k][j].setStress(3,cells[k][j].getDev()[3] + cells[k][j].getPressure());
                }
            }
        }

        /**
         * Check if V at img and hard surfaces is not 0
         */

        for(int i = 0; i<y_count; i++){
            if(nodes[1][i].getvV()[0] != 0.0){
                throw new WrongSpeedException("Speed: " + Arrays.toString(nodes[1][i].getvV()) + " Num: " + Arrays.toString(nodes[1][i].getNum()));
            }
        }

        for(int i = 0; i<x_count; i++){
            if(nodes[i][1].getvV()[1] != 0.0){
                throw new WrongSpeedException("Speed: " + Arrays.toString(nodes[i][1].getvV()) + " Num: " + Arrays.toString(nodes[i][1].getNum()));
            }
        }

//        y_dot_next[j][k]=y_dot_curr[j][k] + (dt/(2*phi[j][k]))*(sigma_yy[j][k]*(x[j+1][k]-x[j][k+1]) +
//                sigma_yy[j-1][k]*(x[j][k+1]-x[j-1][k])+
//                sigma_yy[j-1][k-1]*(x[j-1][k]-x[j][k-1])+
//                sigma_yy[j][k-1]*(x[j][k-1]-x[j+1][k])-
//
//                sigma_xy[j][k]*(y[j+1][k]-y[j][k+1])-
//                sigma_xy[j-1][k]*(y[j][k+1]-y[j-1][k])-
//                sigma_xy[j-1][k-1]*(y[j-1][k]-y[j][k-1])-
//                sigma_xy[j][k-1]*(y[j][k-1]-y[j+1][k]))+dt*betta[j][k];

//        x_dot_next[j][k]=x_dot_curr[j][k] - (dt/(2*phi[j][k]))*(sigma_xx[j][k]*(y[j+1][k]-y[j][k+1]) +
//                sigma_xx[j-1][k]*(y[j][k+1]-y[j-1][k])+
//                sigma_xx[j-1][k-1]*(y[j-1][k]-y[j][k-1])+
//                sigma_xx[j][k-1]*(y[j][k-1]-y[j+1][k])-
//                sigma_xy[j][k]*(x[j+1][k]-x[j][k+1])-
//                sigma_xy[j-1][k]*(x[j][k+1]-x[j-1][k])-
//                sigma_xy[j-1][k-1]*(x[j-1][k]-x[j][k-1])-
//                sigma_xy[j][k-1]*(x[j][k-1]-x[j+1][k]))+dt*alpha[j][k];

//        betta[j][k] = 0.25*((sigma_yy[j][k]- sigma_zz[j][k])*A[j][k]/mass[j][k] +
//                (sigma_yy[j-1][k]- sigma_zz[j-1][k])*A[j-1][k]/mass[j-1][k] +
//                (sigma_yy[j-1][k-1]- sigma_zz[j-1][k-1])*A[j-1][k-1]/mass[j-1][k-1] +
//                (sigma_yy[j][k-1]- sigma_zz[j][k-1])*A[j][k-1]/mass[j][k-1]);

//        alpha[j][k] = 0.25*(sigma_xy[j][k]*A[j][k]/mass[j][k] +
//                sigma_xy[j-1][k]*A[j-1][k]/mass[j-1][k] +
//                sigma_xy[j-1][k-1]*A[j-1][k-1]/mass[j-1][k-1] +
//                sigma_xy[j][k-1]*A[j][k-1]/mass[j][k-1]);

//        phi[j][k]=0.25*(A[j][k]*density_0[j][k]/v_curr[j][k] +
//                A[j-1][k]*density_0[j-1][k]/v_curr[j-1][k] +
//                A[j-1][k-1]*density_0[j-1][k-1]/v_curr[j-1][k-1]+
//                A[j][k-1]*density_0[j][k-1]/v_curr[j][k-1]);

//        /**
//         * =*=*=*=*=*=*=_ Усредняем скорости _=*=*=*=*=*=*=
//         * TODO: Учесть шаг во времени, граничные условия
//         */
//        System.out.println(ANSI_RED+"---->> Bluring veclocity field"+ANSI_RESET);
//
//        for(int i = 1; i < x_count-1; i++){
//            for(int j = 1; j < y_count-1; j++){
//
//                double[][] mm = maxAndMin(nodes[i][j], true);
//
//                double newVx = nodes[i][j].getV()[0]*(1-viscosity_coff) + viscosity_coff*0.5*(mm[0][0] + mm[1][0]);
//                double newVy = nodes[i][j].getV()[1]*(1-viscosity_coff) + viscosity_coff*0.5*(mm[0][1] + mm[1][1]);
//                nodes[i][j].setV(newVx, newVy);
//            }
//        }
//        System.out.println(ANSI_GREEN+"---->> OK"+ANSI_RESET);
        printFields(true, true, false, true, false);
        writeToFile();
    }

    private double op1(int switcher, int numX, int numY) throws WrongAccesException {
        double y3y1, y2y4, x3x1,x2x4;
        switch (switcher){
            case 0: //Velocity x||x
                y3y1 = (oldNodes[numX + 1][numY+1].getCoords()[1] - oldNodes[numX][numY].getCoords()[1] + nodes[numX + 1][numY+1].getCoords()[1] - nodes[numX][numY].getCoords()[1])*0.5;
                y2y4 = (oldNodes[numX + 1][numY].getCoords()[1] - oldNodes[numX][numY+1].getCoords()[1] + nodes[numX + 1][numY].getCoords()[1] - nodes[numX][numY+1].getCoords()[1])*0.5;
                return ((nodes[numX + 1][numY].getV()[0]-nodes[numX][numY+1].getV()[0])*y3y1 -
                (nodes[numX + 1][numY+1].getV()[0]-nodes[numX][numY].getV()[0])*y2y4);
            case 1: //Velocity y||x
                y3y1 = (oldNodes[numX + 1][numY+1].getCoords()[1] - oldNodes[numX][numY].getCoords()[1] + nodes[numX + 1][numY+1].getCoords()[1] - nodes[numX][numY].getCoords()[1])*0.5;
                y2y4 = (oldNodes[numX + 1][numY].getCoords()[1] - oldNodes[numX][numY+1].getCoords()[1] + nodes[numX + 1][numY].getCoords()[1] - nodes[numX][numY+1].getCoords()[1])*0.5;
                return ((nodes[numX + 1][numY].getV()[1]-nodes[numX][numY+1].getV()[1])*y3y1 -
                        (nodes[numX + 1][numY+1].getV()[1]-nodes[numX][numY].getV()[1])*y2y4);
            case 2: // Velocity x||y
                x3x1 = (oldNodes[numX + 1][numY+1].getCoords()[0] - oldNodes[numX][numY].getCoords()[0] + nodes[numX + 1][numY+1].getCoords()[0] - nodes[numX][numY].getCoords()[0])*0.5;
                x2x4 = (oldNodes[numX + 1][numY].getCoords()[0] - oldNodes[numX][numY+1].getCoords()[0] + nodes[numX + 1][numY].getCoords()[0] - nodes[numX][numY+1].getCoords()[0])*0.5;
                return -((nodes[numX + 1][numY].getV()[0]-nodes[numX][numY+1].getV()[0])*x3x1 -
                        (nodes[numX + 1][numY+1].getV()[0]-nodes[numX][numY].getV()[0])*x2x4);
            case 3: // Velocity y||y
                x3x1 = (oldNodes[numX + 1][numY+1].getCoords()[0] - oldNodes[numX][numY].getCoords()[0] + nodes[numX + 1][numY+1].getCoords()[0] - nodes[numX][numY].getCoords()[0])*0.5;
                x2x4 = (oldNodes[numX + 1][numY].getCoords()[0] - oldNodes[numX][numY+1].getCoords()[0] + nodes[numX + 1][numY].getCoords()[0] - nodes[numX][numY+1].getCoords()[0])*0.5;
                return -((nodes[numX + 1][numY].getV()[1]-nodes[numX][numY+1].getV()[1])*x3x1 -
                        (nodes[numX + 1][numY+1].getV()[1]-nodes[numX][numY].getV()[1])*x2x4);
            default:
                throw new WrongAccesException(String.valueOf(switcher));
        }
    }

    private void copyTimeStep() throws WrongAccesException {
        for(Node[] a: nodes){
            for(Node b:a){
                oldNodes[b.getNum()[0]][b.getNum()[1]].setV(b.getvV());
            }
        }
    }

    private void setAreaAndMass(int k, int j, double mass) throws WrongAccesException {
        if(k<x_count-1 && j<y_count-1) {
            double A_a = (nodes[k + 1][j].getCoords()[0] * (nodes[k + 1][j + 1].getCoords()[1] - nodes[k][j + 1].getCoords()[1]) +
                    nodes[k + 1][j + 1].getCoords()[0] * (nodes[k][j + 1].getCoords()[1] - nodes[k + 1][j].getCoords()[1]) +
                    nodes[k][j + 1].getCoords()[0] * (nodes[k + 1][j].getCoords()[1] - nodes[k + 1][j + 1].getCoords()[1])) * 0.5;

            double A_b = (nodes[k + 1][j].getCoords()[0] * (nodes[k][j + 1].getCoords()[1] - nodes[k][j].getCoords()[1]) +
                    nodes[k][j + 1].getCoords()[0] * (nodes[k][j].getCoords()[1] - nodes[k + 1][j].getCoords()[1]) +
                    nodes[k][j].getCoords()[0] * (nodes[k + 1][j].getCoords()[1] - nodes[k][j + 1].getCoords()[1])) * 0.5;

            double V_tmp = (nodes[k + 1][j].getCoords()[1] + nodes[k + 1][j + 1].getCoords()[1] + nodes[k][j + 1].getCoords()[1]) *
                    A_a + (nodes[k + 1][j].getCoords()[1] + nodes[k][j].getCoords()[1] + nodes[k][j + 1].getCoords()[1]) *
                    A_b;

            cells[k][j].setArea(A_a+A_b);
            oldCells[k][j].setArea(A_a+A_b);

            mass *= V_tmp;
        } else {
            mass = 0;
        }

        cells[k][j].setMass(mass);
        oldCells[k][j].setMass(mass);
    }

    private double[][] maxAndMin(Node node, boolean prevStep) throws WrongAccesException {
        int[] num = node.getNum();
        Node curr;
        Node[][] cnodes;

        if(prevStep) {
            cnodes = oldNodes;
        }else {
            cnodes = nodes;
        }

        /**
         *            # -- 0 +1
         *  -1 0 -- # @ # -- +1 0
         *            # -- 0 -1
         */

        curr = cnodes[num[0]][num[1]+1];
        double minX = curr.getV()[0];
        double maxX = minX;
        double minY = curr.getV()[1];
        double maxY = minY;

        curr = cnodes[num[0]-1][num[1]];
        if(minX > curr.getV()[0]){
            minX = curr.getV()[0];
        }
        if(maxX < curr.getV()[0]){
            maxX = curr.getV()[0];
        }
        if(minY > curr.getV()[1]){
            minY = curr.getV()[1];
        }
        if(maxY < curr.getV()[1]){
            maxY = curr.getV()[1];
        }

        curr = cnodes[num[0]+1][num[1]];
        if(minX > curr.getV()[0]){
            minX = curr.getV()[0];
        }
        if(maxX < curr.getV()[0]){
            maxX = curr.getV()[0];
        }
        if(minY > curr.getV()[1]){
            minY = curr.getV()[1];
        }
        if(maxY < curr.getV()[1]){
            maxY = curr.getV()[1];
        }

        curr = cnodes[num[0]][num[1]-1];
        if(minX > curr.getV()[0]){
            minX = curr.getV()[0];
        }
        if(maxX < curr.getV()[0]){
            maxX = curr.getV()[0];
        }
        if(minY > curr.getV()[1]){
            minY = curr.getV()[1];
        }
        if(maxY < curr.getV()[1]){
            maxY = curr.getV()[1];
        }

        return new double[][]{{maxX, maxY},{minX, minY}};
    }

    private void printFields(boolean vel, boolean coords, boolean flags, boolean cellFlags, boolean mass) throws WrongAccesException {
        if(vel) {
            System.out.println(ANSI_YELLOW + "#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#" + ANSI_RESET);
            System.out.println(ANSI_YELLOW + "#-#-#-#-#-#-#-#-# Velocity Field: #-#-#-#-#-#-#-#-#" + ANSI_RESET);
            System.out.println(ANSI_YELLOW + "#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#" + ANSI_RESET);

            for (int i = y_count - 1; i >= 0; i--) {
                StringBuilder pr = new StringBuilder();
                for (int j = 0; j < x_count; j++) {

                    pr.append(nodes[j][i].getColor());
                    pr.append(nodes[j][i].getBgcolor());

                    pr.append("{").append(nodes[j][i].getvV()[0] < 0 ? "" : "+");
                    pr.append(String.format("%.5f", nodes[j][i].getvV()[0]));
                    pr.append(" || ");
                    pr.append(nodes[j][i].getvV()[1] < 0 ? "" : "+");
                    pr.append(String.format("%.5f", nodes[j][i].getvV()[1]));
                    pr.append("}");

//                    if (j != 0 && i != 0) {
                        pr.append(ANSI_RESET);
//                    }

                    pr.append("  ");
                }
                System.out.println(pr.toString());
            }
        }
        System.out.println("\n\n");
        if(coords) {
            System.out.println(ANSI_GREEN + "#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#" + ANSI_RESET);
            System.out.println(ANSI_GREEN + "#-#-#-#-#-#-#-#-# Coords Field: #-#-#-#-#-#-#-#-#" + ANSI_RESET);
            System.out.println(ANSI_GREEN + "#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#" + ANSI_RESET);

            for (int i = y_count - 1; i >= 0; i--) {
                StringBuilder pr = new StringBuilder();
                for (int j = 0; j < x_count; j++) {

                    pr.append(nodes[j][i].getColor());
                    pr.append(nodes[j][i].getBgcolor());

                    pr.append("{").append(nodes[j][i].getvCoords()[0] < 0 ? "" : "+");
                    pr.append(String.format("%.5f", nodes[j][i].getvCoords()[0]));
                    pr.append(" || ");
                    pr.append(nodes[j][i].getvCoords()[1] < 0 ? "" : "+");
                    pr.append(String.format("%.5f", nodes[j][i].getvCoords()[1]));
                    pr.append("}");

                        pr.append(ANSI_RESET);

                    pr.append("  ");
                }
                System.out.println(pr.toString());
            }
        }

        if(flags) {
            System.out.println(ANSI_GREEN + "#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#" + ANSI_RESET);
            System.out.println(ANSI_GREEN + "#-#-#-#-#-#-#-#-# Flags Field: #-#-#-#-#-#-#-#-#" + ANSI_RESET);
            System.out.println(ANSI_GREEN + "#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#" + ANSI_RESET);

            for (int i = y_count - 1; i >= 0; i--) {
                StringBuilder pr = new StringBuilder();
                for (int j = 0; j < x_count; j++) {

                    pr.append(nodes[j][i].getBgcolor()).append(ANSI_BLACK);

                    pr.append("{");
                    pr.append(nodes[j][i].getFlag());
                    pr.append("}");

                    pr.append(ANSI_RESET);

                    pr.append("  ");
                }
                System.out.println(pr.toString());
            }
        }

        if(cellFlags) {
            System.out.println(ANSI_GREEN + "#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#" + ANSI_RESET);
            System.out.println(ANSI_GREEN + "#-#-#-#-#-#-#-#-# Cell Flags Field: #-#-#-#-#-#-#-#-#" + ANSI_RESET);
            System.out.println(ANSI_GREEN + "#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#" + ANSI_RESET);

            for (int i = y_count - 1; i >= 0; i--) {
                StringBuilder pr = new StringBuilder();
                for (int j = 0; j < x_count; j++) {

                    if (cells[j][i].getFlag() == 0) {
                        pr.append(ANSI_POS_BACKGROUND);
                    }

                    pr.append("{");
                    pr.append(cells[j][i].getFlag());
                    pr.append("}");

                    pr.append(ANSI_RESET);

                    pr.append("  ");
                }
                System.out.println(pr.toString());
            }
        }

        if(mass) {
            System.out.println(ANSI_GREEN + "#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#" + ANSI_RESET);
            System.out.println(ANSI_GREEN + "#-#-#-#-#-#-#-#-# Mass Field: #-#-#-#-#-#-#-#-#" + ANSI_RESET);
            System.out.println(ANSI_GREEN + "#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#" + ANSI_RESET);

            for (int i = y_count - 1; i >= 0; i--) {
                StringBuilder pr = new StringBuilder();
                for (int j = 0; j < x_count; j++) {

                    if (cells[j][i].getFlag() == 0) {
                        pr.append(ANSI_POS_BACKGROUND);
                    }

                    pr.append("{");
                    pr.append(String.format("%.5f",cells[j][i].getMass()));
                    pr.append("}");

                    pr.append(ANSI_RESET);

                    pr.append("  ");
                }
                System.out.println(pr.toString());
            }
        }
    }

    private void checkCfl(){

    }

    private static class Cell {
        private double mass;
        private final int[] num;
        private double density = rho;
        private double area = 0;
        /**
         *  flags: 0 - принадлежит телу, 1 - мнимая, 2 - отраженная, 3 4 5 6 - углы, 7 8 9 10 - стенки
         *
         *  1 4 9 9 9 9 9 9 9 9 9 5
         *  1 8 0 0 0 0 0 0 0 0 0 10
         *  1 8 0 0 0 0 0 0 0 0 0 10
         *  1 8 0 0 0 0 0 0 0 0 0 10
         *  1 8 0 0 0 0 0 0 0 0 0 10
         *  1 3 7 7 7 7 7 7 7 7 7 6
         *  1 2 2 2 2 2 2 2 2 2 2 2
         *
         *  stress [xx, yy, xy, tt]
         *  dev [xx, yy, zz]
         */
        private double[] stress = new double[] {0,0,0,0};
        private int flag = 0;
        private double[] dev_stress = new double[] {0,0,0,0};


        Cell(int[] num) {
            this.mass = 0;
            this.num = num;
        }

        Cell(int[] num, double mass) {
            this.mass = mass;
            this.num = num;
        }

        void setMass(double mass) {
            this.mass = mass;
        }

        double getMass() {
            if(flag == 1){
                return oldCells[num[0]+1][num[1]].getMass();
            } else if(flag == 2){
                return oldCells[num[0]][num[1]+1].getMass();
            } else {
                return mass;
            }
        }

        int getFlag() {
            return flag;
        }

        void setFlag(int flag) {
            this.flag = flag;
            if(flag == 1){
                this.density = 0;
                this.mass = 0;
            } else if(flag == 3){
                this.area = 0;
                this.stress = new double[] {0,0,0,0};
                this.area = 0;
                this.area = 0;
            }
        }

        public double getPressure(){
            return (dev_stress[0]+ dev_stress[1]+ dev_stress[2])/3.0;
        }

        public void setDensity(double density) {
            if(this.flag != 1) {
                this.density = density;
            }
        }

        public double getDensity() {
            return density;
        }

        public void setDev(int component, double val) {
            this.dev_stress[component] = val;
        }

        public double[] getDev() {
            return dev_stress;
        }

        public void setStress(int component, double val) {
            this.stress[component] = val;
        }

        public double[] getStress(){
            if(num[1]==0){
                double[] t = oldCells[num[0]][num[1]+1].getStress();
                t[1] *= -1;
                t[2] *= -1;
                t[3] = stress[3];
                if(t[1] != 0){
                    System.out.println("!! ACHTUNG !! Non-Zero stress yy value at cell " + num[0]);
                }
                if(t[3] != 0){
                    System.out.println("!! ACHTUNG !! Non-Zero stress tt value at cell " + num[0]);
                }
                return t;
            } else if(num[0] == 0){
                double[] t = oldCells[num[0]+1][num[1]].getStress();
                t[2] *= -1;
                return t;
            } else if(num[1] == y_count-1 || num[0] == x_count -1){
                return new double[]{0,0,0,0};
            } else {
                return stress;
            }
        }

        public int[] getNum() {
            return num;
        }

        public double getArea() {
//            if()
            if(flag == 1){
                return oldCells[num[0]+1][num[1]].getArea();
            } else {
                return area;
            }
        }

        public void setArea(double area) {
            this.area = area;
        }

        public double getRelDensity(){
            return rho/this.density;
        }
    }

    private static class Node {
        private String color = ANSI_BLACK;;
        private String bgcolor = ANSI_CYAN_BACKGROUND;
        private double[] coords;
        private double[] v;
        private final int[] num;
        private int flag = 0;
//        double phi;
//        double alpha;
//        double beta;

        Node(double[] coords, double[] v, int[] num) {
            this.coords = coords;
            this.v = v;
            this.num = num;
//            this.phi = phi;
//            this.alpha = alpha;
//            this.beta = beta;
        }

        void setV(double vx, double vy) {
            this.v = new double[] {vx, vy};
        }

        void setV(double[] xy) {
            this.v = xy;
        }

        void setCoords(double x, double y) {
            this.coords = new double[] {x, y};
        }

        double[] getV() throws WrongAccesException {
            if(flag == 1 || flag == 2) {
                throw new WrongAccesException();
            }
//            } else if(num[0] == 1){
//                return new double[]{0,0};
//            }
            return v;
        }

        double[] getvV() {
//            if(num[0] == 1){
//                return new double[]{0,0};
//            }
            return v;
        }


        int[] getNum(){
            return num;
        }

        double[] getCoords() throws WrongAccesException{
            if(num[1] == 0 && num[0] == 0){
                throw new WrongAccesException();
            }

            if(num[1] == 0){
                return new double[]{oldNodes[num[0]][num[1]+2].getCoords()[0], oldNodes[num[0]][num[1]+2].getCoords()[1]*-1.0};
            } else if(flag == 1){
                return new double[]{oldNodes[num[0]+2][num[1]].getCoords()[0]*-1.0, oldNodes[num[0]+2][num[1]].getCoords()[1]};
            }else {
                return coords;
            }
        }

        double[] getvCoords(){
            if(num[1] == 0){
                return new double[]{nodes[num[0]][num[1]+2].getvCoords()[0], nodes[num[0]][num[1]+2].getvCoords()[1]*-1.0};
            } else if(flag == 1){
                return new double[]{nodes[num[0]+2][num[1]].getvCoords()[0]*-1.0, nodes[num[0]+2][num[1]].getvCoords()[1]};
            }else {
                return coords;
            }
        }

        int getFlag() {
            return flag;
        }

        void setFlag(int flag) {
            this.flag = flag;
            if(flag == 1){
                this.color = ANSI_BLACK;
                this.bgcolor = ANSI_WHITE_BACKGROUND;
            } else if(flag == 8){
                this.color = ANSI_BLACK;
                this.bgcolor = ANSI_YELLOW_BACKGROUND;
            } else if(flag == 9){
                this.color = ANSI_BLACK;
                this.bgcolor = ANSI_GREEN_BACKGROUND;
            } else if(flag == 10){
                this.color = ANSI_BLACK;
                this.bgcolor = ANSI_PURPLE_BACKGROUND;
            } else if(flag == 7){
                this.color = ANSI_BLACK;
                this.bgcolor = ANSI_BLUE_BACKGROUND;
            } else if(flag == 2){
                this.color = ANSI_BLACK;
                this.bgcolor = ANSI_RED_BACKGROUND;
            }
        }

        public void setColor(String color) {
            this.color = color;
        }

        public String getColor() {
            return color;
        }

        public void setBgcolor(String bgcolor) {
            this.bgcolor = bgcolor;
        }

        public String getBgcolor() {
            return bgcolor+"\u001b[1m";
        }
    }

    private static class WrongAccesException extends Throwable{
        public WrongAccesException() { super(); }
        public WrongAccesException(String message) { super(message); }
        public WrongAccesException(String message, Throwable cause) { super(message, cause); }
        public WrongAccesException(Throwable cause) { super(cause); }
    }

    private static class WrongSpeedException extends Throwable{
        public WrongSpeedException() { super(); }
        public WrongSpeedException(String message) { super(message); }
        public WrongSpeedException(String message, Throwable cause) { super(message, cause); }
        public WrongSpeedException(Throwable cause) { super(cause); }
    }

    private class ImgNode extends Node {
        ImgNode(double[] coords, double[] v, int[] num) {
            super(coords, v, num);
        }
    }

}
