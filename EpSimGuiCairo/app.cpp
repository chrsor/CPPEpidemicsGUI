/**
 * Epidemic-Simulation-Gui
 * @author Christopher Sorg
 * In collaboration with @author Vjekoslav Drvar
 * who implemented the parts with mark 'Simulation_Generator'
 *
 * The program gets NxN matrices
 * m is the length of one side of a partition (by loading a file)
 * (so one m x m partition will be generated to one pixel by Cairomm)
 * One matrix stands for one day in the simulation
*/

//Commented includes are for (future) Windows version
//random() needs to be changed to rand() then
//#include <windows.h>
//#include <shellapi.h>
//#include <algorithm>
#include <gtkmm.h>
#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <optional>
#include <cstdlib>
#include <list>
#include <cmath>
#include <ctime>
#include <set>

using namespace std;

/**
 * Functions from Simulation_Generator
 * @author Vjekoslav Drvar
 */

//random sign generator
int rand_sign() {
    auto r = double(random() % 2);
    return int((r * 2 - 1));
}

//random alpha distance generator
double rand_alpha(double alpha) {
    if (alpha <= 1) return -1;

    auto r = (double) random() / RAND_MAX;
    return pow((1 - r), 1 / (1 - alpha));
}

//random exponential distributed number generator
double rand_exp(double lambda) {
    if (lambda <= 0) return -1;

    auto r = (double) random() / RAND_MAX;
    return -log(1 - r) / lambda;
}

//random poisson distributed number generator
int rand_poisson(const std::vector<double> &v) {
    auto r = (double) random() / RAND_MAX;
    return int(lower_bound(v.begin(), v.end(), r) - v.begin() - 1);
}


//random horizontal move in the grid
int rand_a(int D) {
    return int(random() % (D + 1));
}

void simulate(int N_sim,double Alpha,double Lambda_Poi,double Lambda_Exp_inc,double Lambda_Exp_rec,int Days){
    typedef std::pair<int, int> pairs;

    std::ofstream outfile("TESTsim.txt");

    std::vector<int> aux;
    aux.resize(2);

    std::string yn;

    int n_d = N_sim;
    double dist_alpha_d = Alpha;
    double lambda_poisson_d = Lambda_Poi;
    double lambda_exp_inc_d = Lambda_Exp_inc;
    double lambda_exp_rec_d = Lambda_Exp_rec;
    int days_d = Days;

    //poisson specs
    const int n_poisson = 33;
    vector<double> cdf_poisson;

    cdf_poisson.resize(n_poisson + 1);

    cdf_poisson[0] = 0;
    cdf_poisson[1] = exp(-lambda_poisson_d);
    cdf_poisson[n_poisson] = 1;

    for (int i = 2; i < n_poisson; i++) {
        cdf_poisson[i] = cdf_poisson[i - 1] * lambda_poisson_d / (i - 1);
    }

    double cumulator = 0;

    for (int i = 0; i < cdf_poisson.size() - 1; i++) {
        cumulator += cdf_poisson[i];
        cdf_poisson[i] = cumulator;
    }

    //information matrices/vectors
    vector<vector<int> > grid;
    vector<vector<int> > incubation;
    vector<vector<int> > recovery;
    set<pairs> exposed;
    set<pairs> infectious;
    set<pairs> rem_exposed;
    set<pairs> rem_infectious;

    vector<vector<int> > moves;

    grid.resize(n_d);
    incubation.resize(n_d);
    recovery.resize(n_d);
    moves.resize(0);

    srand(time(NULL));

    //initialize the grid
    for (int i = 0; i < n_d; i++) {
        grid[i].resize(n_d);
        for (int j = 0; j < n_d; j++) {
            grid[i][j] = 0;
        }
    }

    //initialize the incubation matrix
    for (int i = 0; i < n_d; i++) {
        incubation[i].resize(n_d);
        for (int j = 0; j < n_d; j++) {
            incubation[i][j] = round(rand_exp(lambda_exp_inc_d));
        }
    }

    //initialize the recovery matrix
    for (int i = 0; i < n_d; i++) {
        recovery[i].resize(n_d);
        for (int j = 0; j < n_d; j++) {
            recovery[i][j] = round(rand_exp(lambda_exp_rec_d));
        }
    }

    //initialize the infectious list
    //int patient_zero_x=rand_a(n-1), patient_zero_y=rand_a(n-1);
    int patient_zero_x = (n_d - 1) / 2, patient_zero_y = (n_d - 1) / 2;

    grid[patient_zero_x][patient_zero_y] = 1;
    incubation[patient_zero_x][patient_zero_y] = -1;

    infectious.insert(make_pair(patient_zero_x, patient_zero_y));

    //SIMULATION

    double reproduction, rep_end = 0.0;

    outfile << n_d << endl;

    for (int i = 0; i < n_d; i++) {
        for (int j = 0; j < n_d; j++) {

            char current = 's';

            if (grid[i][j]) {
                if (incubation[i][j]==-1) {
                    current = 'i';
                }else{
                    current='e';
                }
            }

            else if (recovery[i][j]==-1) {
                current = 'r';
            }
            outfile << current << " ";
        }
        outfile << endl;
    }
    outfile << endl;


    for (int day = 1; day <= days_d; day++) {

        reproduction = 0;

        //moves vector setup

        moves.resize(0);

        for (auto const &x : infectious) {

            int l1 = rand_poisson(cdf_poisson);
            int l2 = rand_poisson(cdf_poisson);

            for (int j = 0; j < l1 + l2; j++) {
                aux[0] = x.first;
                aux[1] = x.second;
                moves.push_back(aux);
            }
        }

        random_shuffle(moves.begin(), moves.end());

        int D, a, b, x, y, i2, j2;
        int N = n_d;

        //moves
        for (int z=0;z<moves.size();z++) {

            i2 = moves[z][0];
            j2 = moves[z][1];

            D = round(rand_alpha(dist_alpha_d));

            a = rand_a(D);
            b = D - a;
            x = (i2 + rand_sign() * a) % N;
            y = (j2 + rand_sign() * b) % N;

            if (x < 0) x += N;
            if (y < 0) y += N;

            if ((grid[x][y] == 0) && (recovery[x][y] != -1)) {
                grid[x][y] = 1;
                exposed.insert(make_pair(x, y));
                reproduction += 1;
            }

        }

        for (int i = 0; i < n_d; i++) {
            for (int j = 0; j < n_d; j++) {

                char current = 's';

                if (grid[i][j]) {
                    if (incubation[i][j] == -1) current = 'i';
                    else current = 'e';
                } else if (recovery[i][j] == -1) current = 'r';
                outfile << current << " ";
            }
            outfile << endl;
        }
        outfile << endl;

        for (auto const &x : exposed) {
            incubation[x.first][x.second]-=1;
            if (incubation[x.first][x.second]==-1) {
                rem_exposed.insert(x);
                infectious.insert(x);
            }
        }

        for (auto const &x : rem_exposed) {
            exposed.erase(x);
        }

        rem_exposed.clear();

        for (auto const &x : infectious) {
            recovery[x.first][x.second]-=1;
            if (recovery[x.first][x.second]==-1) {
                rem_infectious.insert(x);
                grid[x.first][x.second]=0;
            }
        }

        if (rem_infectious.size()) reproduction=reproduction/rem_infectious.size(); else reproduction=reproduction;
        rep_end+=reproduction;

        for (auto const &x : rem_infectious) {
            infectious.erase(x);
        }

        rem_infectious.clear();

    }
}
//Functions from Simulation_Generator End

set<char> changedChars;
vector<string> newColors;

struct color_t {
    double red;
    double green;
    double blue;

    /**
     * Coloring chars for the SEIR model
     * For other chars the user will define the color clicking on 'Change Color' keyboarding the RGB-values
     * Attention! RGB-values are normalised from 0 to 1 in Gtk::Cairo instead from 0 to 255 (default)
     * @param c
     * @return color_t of char c
     */

    static color_t from_char(char c) {
        switch (c) { //RGB Colors, so color_t with {(double)red,(double)green,(double)blue}
            case 's':
                return color_t{0, 0.5, 0}; //Color green
            case 'e':
                return color_t{1, 0.5, 0}; //Color orange
            case 'i':
                return color_t{1, 0, 0}; //Color red
            case 'r':
                return color_t{0, 0, 1}; //Color blue
            case 'd':
            default:
                return color_t{0, 0, 0}; //Color black
        }
    }
};

/**
 * Overloading the + Operator for colors by adding the RGB-values itself
 * @param base , first summand
 * @param calc , second summand
 * @return color_t
 */

color_t operator+(const color_t &base, const color_t &calc) {
    return color_t{
            base.red + calc.red,
            base.green + calc.green,
            base.blue + calc.blue
    };
}

/**
 * Overloading the / Operator for colors by dividing the RGB-values itself
 * @param base , dividend
 * @param i , divisor
 * @return color_t
 */

color_t operator/(const color_t &base, int i) {
    return color_t{
            base.red / i,
            base.green / i,
            base.blue / i
    };
}

struct bitmap_t {
    std::vector<char> data;
    int N{}; //First line of txt file

    /**
     * In dependence of existence of changed char-colors this method
     * either will get the saved RGB-values the user has keyboarded
     * (those chars are saved in set changedChars, the char with its color values is saved in the list newColors)
     * or will color the chars with the defined colors in method from_char
     * (All chars which are not element of {s,e,i,r} will be colored black by default)
     * @param x , cartesian coordinate x
     * @param y , cartesian coordinate y
     * @param m , chosen partition length
     * @return color_t color_mean, the mean value of the mxm partition
     */

    color_t get_color(int x, int y, int m) {

        m = std::min(m, N);
        auto color_sum = color_t::from_char('d'); //To get the default value (0,0,0)

        for (int i = y; i < m + y; i++) {
            for (int j = x; j < m + x; j++) {
                char current = data[i * N + j];
                if (changedChars.contains(current)) {
                    for (int k = 0; k + 4 <= newColors.size(); k += 4) {
                        if (newColors[k].at(0) == current) {

                            double r = atof(newColors.at(k + 1).c_str());
                            double g = atof(newColors.at(k + 2).c_str());
                            double b = atof(newColors.at(k + 3).c_str());

                            auto color = color_t{r / 255, g / 255, b / 255};
                            color_sum = color_sum + color;
                            break;
                        }
                    }
                } else {
                    auto color = color_t::from_char(current);
                    color_sum = color_sum + color; //Defined + Operator
                }
            }
        }
        auto color_mean = color_sum / (m * m); //Defined / Operator
        return color_mean;
    }

    /**
     * Reading the file and saving it in vector data
     * @param input
     * @param N
     * @return bmp
     */

    static std::optional<bitmap_t> from_stream(std::istream &input, int N) {
        bitmap_t bmp;
        bmp.data.reserve(N * N); //Reserve N*N space in bmp (similar to functions like malloc)
        bmp.N = N;
        char c;

        for (int i = 0; i < (N) * (N); i++) {
            if (input.eof()) { //Reading-Error
                return {};
            }
            input >> c;
            if (c != '\0') {
                bmp.data.push_back(c);
            }
        }
        return bmp;
    }

};

/**
 * Validating data and loading
 * @param input
 * @return bit_valid
 */

std::vector<bitmap_t> load(std::istream &input) {
    int N;
    input >> N;
    if (N <= 0) { //A matrix M(NxN,R) with N <= 0 doesn't make any sense (here)
        return {};
    }

    std::vector<bitmap_t> bit_valid;
    std::optional<bitmap_t> bit;
    while ((bit = bitmap_t::from_stream(input, N))) {
        bit_valid.push_back(bit.value());
    }

    return bit_valid;
}


std::vector<bitmap_t> graphical_bitmap_vector;
int graphical_m;
int graphical_day;
double Scaling_Factor = 10.00;

/**
 * Reseting values for m and day and clearing the graphical_bitmap_vector
 * So the default values are: m=1, day=0
 */

void graphical_reset() {
    graphical_bitmap_vector.clear();
    graphical_m = 1;
    graphical_day = 0;
}

/**
 * Initialise the slider for the GUI
 * @param slider
 */

void initialise_slider(Gtk::Scale *slider){
    Glib::RefPtr<Gtk::Adjustment> g_adjustment;
    g_adjustment = Gtk::Adjustment::create(0, 0, double(graphical_bitmap_vector.size() - 1), 1);
    slider->set_adjustment(g_adjustment);
    graphical_day = int(slider->get_value());
}

/**
 * Take screenshot of the current state of the GUI
 * (for zooming in the picture of special day for example)
 */

void take_screenshot(){
    GdkWindow *source;
    GdkPixbuf *screenshot;
    gint x, y, width, height;

    source = gdk_get_default_root_window();
    gdk_window_get_geometry (source, &x, &y, &width, &height);
    screenshot = gdk_pixbuf_get_from_window (source, x, y, width, height);

    gdk_pixbuf_save(screenshot,"chosenDay.png","png",nullptr,NULL);
}

/**
 * Method for dealing with all buttons and widgets in the GUI
 * The method uses lambda capture lists ([=] and []) for connecting the buttons to the right slot
 * @param app
 */

void on_app_activate(const Glib::RefPtr<Gtk::Application> &app) {
    Gtk::Window *window = nullptr;
    Gtk::DrawingArea *drawar = nullptr;
    Gtk::Button *info_button = nullptr;
    Gtk::Button *load_button = nullptr;
    Gtk::Button *scale_button = nullptr;
    Gtk::Button *allColor_button = nullptr;
    Gtk::Button *m_button = nullptr;
    Gtk::Button *day_button = nullptr;
    Gtk::Button *sim_button = nullptr;
    Gtk::Button *rgb_button = nullptr;
    Gtk::Button *color_button = nullptr;
    Gtk::Button *screenshot_button = nullptr;
    Gtk::Scale *dayscale = nullptr;
    Gtk::Label *info_label = nullptr;


    auto builder = Gtk::Builder::create_from_file("../window.glade");
    builder->get_widget("window", window);
    builder->get_widget("drawar", drawar);
    builder->get_widget("info_button", info_button);
    builder->get_widget("load_button", load_button);
    builder->get_widget("scale_button", scale_button);
    builder->get_widget("m_button", m_button);
    builder->get_widget("rgb_button", rgb_button);
    builder->get_widget("day_button", day_button);
    builder->get_widget("sim_button", sim_button);
    builder->get_widget("info_label", info_label);
    builder->get_widget("color_button", color_button);
    builder->get_widget("screenshot_button",screenshot_button);
    builder->get_widget("dayscale", dayscale);
    builder->get_widget("allColor_button", allColor_button);

    //Slot [=] is a lambda capture list with Copy-By-Value
    auto update_info = [=]() {
        info_label->set_text(
                "PartitionLength: " + std::to_string(graphical_m) + " || Scalar: " + std::to_string(Scaling_Factor) +
                " || Day: " + std::to_string(graphical_day));
    };

    //Message that occurs when clicking on 'Help'
    info_button->signal_clicked().connect([]() {
        Gtk::MessageDialog("Hello! This is an Epidemic Simulation GUI.\n\n"
                           "Click on 'Load File' to choose a txt file in SEIR-Matrix-Format to load.\n"
                           "(An example for the SEIR-Matrix-Format can be found by opening the file 'TEST.txt')\n\n"
                           "Click on 'Change Partition' to change the length of the m x m Partitions.\n\n"
                           "Click on 'Change Scalar' to change the Scaling Factor of the picture.\n"
                           "Highly recommended: Choose Scalar my, so that 700*700 Px >= (my * N)*(my * N) Px\n\n"
                           "Click on 'Choose Day' to choose the Day to load from the txt file, beginning at Day 0.\n\n"
                           "You can change the day by moving the slider.\n\n"
                           "Click on 'Change Color' to define your own color for a char.\n"
                           "First input is the char, then RGB values (in this order, please use integer values).\n"
                           "If you enter values higher than 255, the RGB values will always be scaled to 255.\n"
                           "If you enter invalid values with format char + numbers, you will define 0 for that value.\n\n"
                           "Click on 'Simulate' to simulate an own Epidemic-Run.\n"
                           "Parameters are: Matrix-Size N, Stochastic Parameters alpha, lambda(Poisson_dist), lambda(Exp_dist_inc), lambda(Exp_dist_rec), Days\n"
                           "ATTENTION! You generate a file called 'Textsim.txt'. If you simulate again, you will overwrite this file if not renamed manually before.\n\n"
                           "Click on 'RGB-Calculator' for opening a website that is able to calculate colors with RGB-values.\n"
                           "ATTENTION! To use this button please open your default browser FIRST and then click on the button.\n"
                           "If not the GUI will only work again if you close the window in your browser.\n\n"
                           "Click on 'Show all colors' to get an overview of the defined colors for all chars.\n"
                           "(You can see there if your own defined color for a char has been saved or not.)\n\n"
                           "Click on 'Screenshot' to save the current state of the GUI as screenshot.").run();
    });

    //Action when clicking on 'Load File'
    load_button->signal_clicked().connect([=]() {
        auto dialog = Gtk::FileChooserDialog(*window, "Choose File");
        dialog.add_button("Select", Gtk::ResponseType::RESPONSE_OK);
        if (dialog.run() != Gtk::ResponseType::RESPONSE_OK) {
            return;
        }

        auto filename = dialog.get_file()->get_path();
        std::ifstream input(filename);
        graphical_reset();
        graphical_bitmap_vector = load(input);

        set<char> inputChars;
        for (int i = 0; i < std::ifstream::end; i++) {
            if (!(inputChars.find(i) != inputChars.end())) {
                inputChars.insert(i);
            }
        }

        initialise_slider(dayscale);
        update_info();
        drawar->queue_draw();
    });

    //Action when clicking on 'Change Partition'
    m_button->signal_clicked().connect([=]() {
        if (graphical_bitmap_vector.empty()) {
            return;
        }
        int N = graphical_bitmap_vector.front().N;

        auto dialog = Gtk::Dialog("Choose m", *window, true);
        Gtk::Entry m_entry;
        dialog.get_content_area()->add(m_entry);
        m_entry.set_text("m");
        m_entry.show();
        dialog.add_button("Ok", Gtk::ResponseType::RESPONSE_OK);
        dialog.add_button("Cancel", Gtk::ResponseType::RESPONSE_CANCEL);
        if (dialog.run() != Gtk::ResponseType::RESPONSE_OK) {
            return;
        }

        auto m_text = m_entry.get_text();
        int m = atoi(m_text.c_str()); //atoi converts String to Integer
        if (m <= 0 || m > N) {
            return;
        }

        graphical_m = m;
        update_info();
        drawar->queue_draw();
    });

    //Action when clicking on 'Change Scalar'
    scale_button->signal_clicked().connect([=]() {
        auto dialog = Gtk::Dialog("Choose Scalar", *window, true);
        Gtk::Entry scale_entry;
        dialog.get_content_area()->add(scale_entry);
        scale_entry.set_text("Scalar");
        scale_entry.show();
        dialog.add_button("Ok", Gtk::ResponseType::RESPONSE_OK);
        dialog.add_button("Cancel", Gtk::ResponseType::RESPONSE_CANCEL);
        if (dialog.run() != Gtk::ResponseType::RESPONSE_OK) {
            return;
        }

        auto scale_text = scale_entry.get_text();
        double scale = atof(scale_text.c_str()); //atof converts String to Double
        if (scale <= 0) {
            return;
        }

        Scaling_Factor = scale;
        update_info();
        drawar->queue_draw();
    });


    //Action when clicking on 'Choose Day'
    day_button->signal_clicked().connect([=]() {
        if (graphical_bitmap_vector.empty()) {
            return;
        }

        auto dialog = Gtk::Dialog("Choose Day", *window, true);

        Gtk::Entry day_entry;
        dialog.get_content_area()->add(day_entry);
        day_entry.set_text("Day");
        day_entry.show();
        dialog.add_button("Ok", Gtk::ResponseType::RESPONSE_OK);
        dialog.add_button("Cancel", Gtk::ResponseType::RESPONSE_CANCEL);
        if (dialog.run() != Gtk::ResponseType::RESPONSE_OK) {
            return;
        }

        auto text = day_entry.get_text();
        int day = atoi(text.c_str());
        if (day < 0 || day >= graphical_bitmap_vector.size()) {
            return;
        }

        graphical_day = day;
        update_info();
        drawar->queue_draw();
    });

    //Action when clicking on 'Simulate'
    sim_button->signal_clicked().connect([=]() {
        auto dialog = Gtk::Dialog("Choose Parameters", *window, true);

        Gtk::Entry N_entry;
        dialog.get_content_area()->add(N_entry);
        N_entry.set_text("N");
        N_entry.show();

        Gtk::Entry Alpha_entry;
        dialog.get_content_area()->add(Alpha_entry);
        Alpha_entry.set_text("Alpha");
        Alpha_entry.show();

        Gtk::Entry Lambda_Poi_entry;
        dialog.get_content_area()->add(Lambda_Poi_entry);
        Lambda_Poi_entry.set_text("Lambda Poisson");
        Lambda_Poi_entry.show();

        Gtk::Entry Lambda_Exp_inc_entry;
        dialog.get_content_area()->add(Lambda_Exp_inc_entry);
        Lambda_Exp_inc_entry.set_text("Lambda Exp Inc");
        Lambda_Exp_inc_entry.show();

        Gtk::Entry Lambda_Exp_rec_entry;
        dialog.get_content_area()->add(Lambda_Exp_rec_entry);
        Lambda_Exp_rec_entry.set_text("Lambda Exp Rec");
        Lambda_Exp_rec_entry.show();

        Gtk::Entry Days_entry;
        dialog.get_content_area()->add(Days_entry);
        Days_entry.set_text("Days");
        Days_entry.show();

        dialog.add_button("Ok", Gtk::ResponseType::RESPONSE_OK);
        dialog.add_button("Cancel", Gtk::ResponseType::RESPONSE_CANCEL);
        if (dialog.run() != Gtk::ResponseType::RESPONSE_OK) {
            return;
        }

        auto N_sim_text = N_entry.get_text();
        int N_sim = atoi(N_sim_text.c_str());

        auto Alpha_sim_text = Alpha_entry.get_text();
        double Alpha_sim = atof(Alpha_sim_text.c_str());

        auto Lambda_Poi_sim_text = Lambda_Poi_entry.get_text();
        double Lambda_Poi_sim = atof(Lambda_Poi_sim_text.c_str());

        auto Lambda_Exp_inc_sim_text = Lambda_Exp_inc_entry.get_text();
        double Lambda_Exp_inc_sim = atof(Lambda_Exp_inc_sim_text.c_str());

        auto Lambda_Exp_rec_sim_text = Lambda_Exp_rec_entry.get_text();
        double Lambda_Exp_rec_sim = atof(Lambda_Exp_rec_sim_text.c_str());

        auto Days_sim_text = Days_entry.get_text();
        int Days_sim = atoi(Days_sim_text.c_str());

        simulate(N_sim,Alpha_sim,Lambda_Poi_sim,Lambda_Exp_inc_sim,Lambda_Exp_rec_sim,Days_sim);

        std::ifstream input("TESTsim.txt");
        graphical_reset();
        update_info();
        graphical_bitmap_vector = load(input);

        set<char> inputChars;
        for (int i = 0; i < std::ifstream::end; i++) {
            if (!(inputChars.find(i) != inputChars.end())) {
                inputChars.insert(i);
            }
        }

        initialise_slider(dayscale);
        drawar->queue_draw();
    });

    //Action when clicking on 'RGB-Calculator'
    rgb_button->signal_clicked().connect([=]() {
        //For Linux
        system("sensible-browser -new-tab https://www.w3schools.com/colors/colors_rgb.asp");
        //For Windows
        //ShellExecute(nullptr, "open", "https://www.w3schools.com/colors/colors_rgb.asp", nullptr, nullptr, SW_SHOW);
    });

    //Action when clicking on 'Change Color'
    color_button->signal_clicked().connect([=]() {
        auto dialog = Gtk::Dialog("Choose Char + RGB", *window, true);

        Gtk::Entry char_entry;
        dialog.get_content_area()->add(char_entry);
        char_entry.set_text("C");
        char_entry.set_max_length(1);
        char_entry.show();

        Gtk::Entry red_entry;
        dialog.get_content_area()->add(red_entry);
        red_entry.set_text("R");
        red_entry.set_max_length(3);
        red_entry.show();

        Gtk::Entry green_entry;
        dialog.get_content_area()->add(green_entry);
        green_entry.set_text("G");
        green_entry.set_max_length(3);
        green_entry.show();

        Gtk::Entry blue_entry;
        dialog.get_content_area()->add(blue_entry);
        blue_entry.set_text("B");
        blue_entry.set_max_length(3);
        blue_entry.show();

        dialog.add_button("Ok", Gtk::ResponseType::RESPONSE_OK);
        dialog.add_button("Cancel", Gtk::ResponseType::RESPONSE_CANCEL);
        if (dialog.run() != Gtk::ResponseType::RESPONSE_OK) {
            return;
        }

        string char_text = char_entry.get_text();
        char c = char_text[0];

        auto red_text = red_entry.get_text();

        auto green_text = green_entry.get_text();

        auto blue_text = blue_entry.get_text();

        //Ascii values for chars begin at 65(A), so small a has Ascii value 97 for example
        //But: You can only tip in one char in char_entry, so if the user is keyboarding a number
        //instead of a char, the value can't be that high
        //At least one of the 3 digits you can enter in RGB should be a number n, so 48 (0) <= n <= 57 (9)
        if ((c >= 65) && (
                (((48 <= red_text[0]) && (red_text[0] <= 57)) || ((48 <= red_text[1]) && (red_text[1] <= 57))
                 || ((48 <= red_text[2]) && (red_text[2] <= 57))) &&
                (((48 <= green_text[0]) && (green_text[0] <= 57)) || ((48 <= green_text[1]) && (green_text[1] <= 57))
                 || ((48 <= green_text[2]) && (green_text[2] <= 57))) &&
                (((48 <= blue_text[0]) && (blue_text[0] <= 57)) || ((48 <= blue_text[1]) && (blue_text[1] <= 57))
                 || ((48 <= blue_text[2]) && (blue_text[2] <= 57)))
        )) {
            //Then the values for that char have to be overwritten
            if (changedChars.contains(c)) {
                for (int k = 0; k + 4 <= newColors.size(); k += 4) {
                    if (newColors[k].at(0) == c) {
                        newColors[k + 1].replace(newColors[k + 1].begin(), newColors[k + 1].end(), red_text);
                        newColors[k + 2].replace(newColors[k + 2].begin(), newColors[k + 2].end(), green_text);
                        newColors[k + 3].replace(newColors[k + 3].begin(), newColors[k + 3].end(), blue_text);
                        break;
                    }
                }
            } else {
                changedChars.insert(c);

                newColors.push_back(char_text);
                newColors.push_back(red_text);
                newColors.push_back(green_text);
                newColors.push_back(blue_text);
            }
        } else {
            auto errorDialog = Gtk::MessageDialog(*window, "Invalid input, char has not been saved.\nPlease try again",
                                                  false,
                                                  Gtk::MessageType::MESSAGE_ERROR,
                                                  Gtk::ButtonsType::BUTTONS_OK, true);
            errorDialog.run();
        }

        update_info();
        drawar->queue_draw();
    });

    //Action when clicking on 'Show all Colors'
    allColor_button->signal_clicked().connect([=]() {
        string msg = "Char Red Green Blue\n\n"
                     "---[SEIR]---\n\n"
                     "s\t0\t128\t0\n"
                     "e\t255\t128\t0\n"
                     "i\t255\t0\t0\n"
                     "r\t0\t0\t255\n\n"
                     "---[NEW]---\n\n";
        if (changedChars.empty()) {
            msg.append("No chars have been defined or changed.");
            auto msgDialog = Gtk::MessageDialog(*window, msg, false,
                                                Gtk::MessageType::MESSAGE_INFO,
                                                Gtk::ButtonsType::BUTTONS_OK, true);
            msgDialog.run();
        } else {
            for (int k = 0; k + 4 <= newColors.size(); k += 4) {
                msg.append(newColors.at(k)).append("\t").append(newColors.at(k + 1)).append("\t").
                        append(newColors.at(k + 2)).append("\t").append(newColors.at(k + 3)).append("\n");
            }
            auto msgDialog = Gtk::MessageDialog(*window, msg, false,
                                                Gtk::MessageType::MESSAGE_INFO,
                                                Gtk::ButtonsType::BUTTONS_OK, true);
            msgDialog.run();
        }
    });

    //Implementing the dayscale slider
    dayscale->signal_value_changed().connect([=]() {
        if (graphical_bitmap_vector.empty()) {}
        else {
            graphical_day = int(dayscale->get_value());
            update_info();
            drawar->queue_draw();
        }
    });

    //Action when clicking on 'Screenshot'
    screenshot_button->signal_clicked().connect([=](){
       take_screenshot();
    });

    //Slot [] means no capturing in lambda
    drawar->signal_draw().connect([](const Cairo::RefPtr<Cairo::Context> &cr) {
        if (graphical_bitmap_vector.empty() || graphical_day < 0 || graphical_m <= 0 ||
            graphical_day >= graphical_bitmap_vector.size()) {
            return false; //False so that Cairo knows to stop drawing
        }
        int N = graphical_bitmap_vector.front().N;
        int parts = N / graphical_m;

        //Coloring of the partitions by Cairo
        for (int i = 0; i < parts; i++) {
            for (int j = 0; j < parts; j++) {
                color_t color = graphical_bitmap_vector[graphical_day].get_color(j * graphical_m, i * graphical_m,
                                                                                 graphical_m);
                cr->save();
                cr->set_source_rgb(color.red, color.green, color.blue);
                cr->rectangle(Scaling_Factor * j, Scaling_Factor * i, Scaling_Factor, Scaling_Factor);
                cr->fill();
                cr->restore();
            }
        }
        return false;
    });

    graphical_reset();
    update_info();
    app->add_window(*window);
    window->show_all();
}

int main(int argc, char **argv) {
    auto app = Gtk::Application::create(argc, argv);

    app->signal_activate().connect([=]() {
        on_app_activate(app);
    });

    return app->run(argc, argv);
}