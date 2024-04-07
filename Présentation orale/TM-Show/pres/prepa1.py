from manim import *
from manim_slides import Slide, ThreeDSlide
import random as rd
import fitz


# manim -qh -a prepa1.py
# manim-slides TitleSlide SommaireSlide RevueInformatiqueSlide PrincipeOrdinateurQuantiqueSlide RevuePhysiqueSlide LanguageMathSlide QubitSlide ExemplesAlgoSlide DeutschSlide GroverSlide SudokuSlide RevueOrdinateurQuantiqueSlide AvenirSlide EndOfPresentation


class TitleSlide(Slide):
    def construct(self):
        title = Text("Les algorithmes quantiques", t2w={"algorithmes quantiques": BOLD})
        subtitle = Text("Ou une théorie d'optimisation", t2s={"optimisation": ITALIC})
        name = Text("Travail de maturité par Romain Blondel", t2c={"Romain": RED, "Blondel": BLUE})
        name_short = Text("Romain Blondel", t2c={"Romain": RED, "Blondel": BLUE})
        supervision = Text("Sous la supervision de M. Frédéric De Montmollin", t2c={"Frédéric": RED, "De Montmollin": BLUE})
        supervision_short = Text("Supervision : M. De Montmollin", t2c={"Frédéric": RED, "De Montmollin": BLUE})
        classe = Text("3M8, Gymnase Auguste Piccard", t2c={"G": YELLOW, "A": YELLOW, "P": YELLOW, "3M8": RED})
        classe_short = Text("3M8, GAP", t2c={"G": YELLOW, "A": YELLOW, "P": YELLOW, "3M8": RED})
        date = Text("5 décembre 2023")
        date_short = Text("5.12.2023")
        title_image = ImageMobject("sd5.png")

        title.scale(1.2)
        subtitle.scale(0.8)
        title_image.scale(1.1)
        title_image.shift(DOWN*0.7)

        for i in [name, supervision, classe, date, name_short, supervision_short, classe_short, date_short]:
            i.scale(0.3)

        title.shift(UP*3)
        subtitle.shift(UP*2)
        name.shift(DOWN*2.5 + LEFT*4)
        name_short.shift(DOWN*2.5 + LEFT*4)
        supervision.shift(DOWN*3 + LEFT*4)
        supervision_short.shift(DOWN*3 + LEFT*4)
        classe.shift(DOWN*3 + RIGHT*4)
        classe_short.shift(DOWN*3 + RIGHT*4)
        date.shift(DOWN*2.5 + RIGHT*4)
        date_short.shift(DOWN*2.5 + RIGHT*4)

        self.play(FadeIn(title, name, supervision, classe, date))

        self.next_slide()

        self.play(Write(subtitle))

        self.next_slide()

        for i in [[name, name_short], [supervision, supervision_short], [classe, classe_short], [date, date_short]]:
            self.play(Transform(i[0], i[1]))

        self.play(GrowFromCenter(title_image))    

        self.next_slide()

        self.play(FadeOut(title, subtitle, title_image, name_short, supervision_short, classe_short, date_short, name, supervision, classe, date))
        

class SommaireSlide(Slide):
    def construct(self):
        title = Text("Sommaire")
        title.shift(UP*3)

        self.play(FadeIn(title))

        self.next_slide()

        sommaire = Tex(r"\begin{itemize} \item[$\bullet$] Revue du travail \\ \begin{itemize} \item Notions de physique et d'informatique \\ \item Language mathématique \\ \item Principe d'un ordinateur quantique \\ \begin{itemize} \item[$\rightarrow$] Parallèle au classique \end{itemize} \item Exemples d'algorithmes \\ \item Perspectives pour l'avenir \end{itemize} \item[$\bullet$] Réalisations concrètes \\ \begin{itemize} \item La physique quantique \\ \item Des qubits \\ \item Algorithme de Deutsch \\ \item Algorithme de Grover - Sudoku \end{itemize}\end{itemize}")
        sommaire.scale(0.5)
        sommaire.shift(LEFT*3 + DOWN*0.5)

        self.play(Write(sommaire))

        self.next_slide()

        self.play(FadeOut(title, sommaire))                


class RevueInformatiqueSlide(Slide):
    def construct(self):
        title = Text("Notions d'informatique")
        title.shift(UP*3)
        title.scale(0.8)

        self.play(FadeIn(title))

        self.next_slide()

        c_code_blank = '''int factorial(int n)
{
    int result = 1;
    for (int i = 1; i <= n; i++)
    {
        result *= i;
    }
    return result;
}
'''

        c_code = Code(code=c_code_blank, language="c", tab_width=4, insert_line_no=False, background="window", style="monokai", font="Monospace")

        python_code_blank = '''def factorial(n):
    result = 1
    for i in range(1, n+1):
        result *= i
    return result
'''

        python_code = Code(code=python_code_blank, language="python", tab_width=4, insert_line_no=False, background="window", style="monokai", font="Monospace")

        c_code.scale(0.6)
        python_code.scale(0.6)

        c_time = Tex(r"Temps d'exécution : ", r"$\sim 10^{-8} \ [s]$")
        c_time.set_color_by_tex("[s]", YELLOW)

        python_time = Tex(r"Temps d'exécution : ", r"$\sim 10^{-5} \ [s]$")
        python_time.set_color_by_tex("[s]", YELLOW)

        c_time.scale(0.8)
        python_time.scale(0.8)

        code = Group(c_code, python_code).arrange(buff=1, direction=DOWN)

        c_time.next_to(c_code, DOWN)
        python_time.next_to(python_code, DOWN)

        self.play(FadeIn(code))

        self.next_slide()

        self.play(FadeIn(c_time))

        self.next_slide()

        self.play(FadeIn(python_time))

        self.next_slide()

        code_example_gpe = Group(code, c_time, python_time)

        self.play(Transform(code_example_gpe, code_example_gpe.shift(LEFT*3)))

        axes_comp = Axes(x_range=[0, 1200000, 1000000], y_range=[0, 350, 100], x_length=4, y_length=5, axis_config={"include_numbers": True})
        axes_comp.shift(RIGHT*3+DOWN*0.5)
        labels = axes_comp.get_axis_labels(
            Tex(r"Nombre\\ en entrée").scale(0.5), Tex(r"Temps de\\ calcul [s]").scale(0.5)
        )

        dots = [[1, 1.59999990501092e-06],
                [100001, 1.80004839999992],
                [200001, 7.41467120000016],
                [300001, 17.3281252000002],
                [400001, 31.7894969000004],
                [500001, 59.6285113999998],
                [600001, 106.5159247],
                [700001, 161.3160416],
                [800001, 226.8987501],
                [900001, 294.3835556]]

        dots_graph = Group()
        for i in dots:
            dots_graph.add(Dot(axes_comp.c2p(i[0], i[1]), color=RED))

        graph_comp = axes_comp.plot(lambda x: (4.2838302e-16 * x**3), x_range=[0, (350 / (4.2838302e-16))**(1/3), 1000], use_smoothing=True, color=YELLOW)

        graph_legend = Tex(r"$a \cdot x^3$")
        graph_legend.next_to(graph_comp, UP+RIGHT, buff=0.1)
        graph_legend.set_color(YELLOW)
        graph_legend.scale(0.7)

        self.play(FadeIn(axes_comp, labels))

        self.next_slide()

        self.play(FadeIn(dots_graph))

        self.next_slide()

        self.play(Write(graph_comp))
        self.play(FadeIn(graph_legend))

        self.next_slide()

        self.play(FadeOut(title, code_example_gpe, dots_graph, graph_comp, graph_legend, axes_comp, labels))


class PrincipeOrdinateurQuantiqueSlide(Slide):
    def construct(self):
        title = Text("La physique quantique")
        title.shift(UP*3)

        self.play(FadeIn(title))

        self.next_slide()

        dec = 3

        light_source = Dot(LEFT*3 + LEFT*dec).set_color(YELLOW)

        wall_with_holes = Group(Line(UP*0.5, DOWN*0.5), Line(UP*0.7, UP*2), Line(DOWN*0.7, DOWN*2))
        wall_with_holes.shift(LEFT*dec)

        screen = Line(DOWN*2, UP*2).shift(RIGHT*3)
        screen.shift(LEFT*dec)

        self.play(FadeIn(light_source, wall_with_holes, screen))

        self.next_slide(auto_next=True)

        self.create_wavefront(-3 + LEFT*dec, 0)

        self.sym_wavefront(0 + LEFT*dec, 0.6)

        luminosity_axes = Axes(x_range=[-20, 20], y_range=[0, 2], x_length=4, y_length=1, axis_config={"include_ticks": False})
        luminosity_axes.rotate(-PI/2)
        luminosity_axes.shift(RIGHT*3 + LEFT*dec)
        luminosity_graph = luminosity_axes.plot(lambda x: np.sin(x)/(x) + 1, x_range=[-20, 20], color=YELLOW)

        self.next_slide(loop=True)

        self.wait(3)

        self.next_slide(auto_next=True)

        self.play(FadeIn(luminosity_graph))

        self.next_slide(loop=True)

        self.wait(3)

        self.next_slide(auto_next=True)

        screen_discreate = RoundedRectangle(height=2, width=4, corner_radius=0.2).shift(RIGHT*3.5, UP*1.5)

        self.play(FadeIn(screen_discreate))
        
        bands = [[-1, 1], [2,3], [3.5, 4], [-2, -3], [-3.5, -4]]
        num_points = [500, 200, 80, 200, 80]

        lights_dot = Group()
        for i in range(len(bands)):
            for j in range(num_points[i]):
                lights_dot.add(Dot(np.random.uniform(bands[i][0]/2.2, bands[i][1]/2.2)*RIGHT + np.random.uniform(-1, 1)*UP + RIGHT*3.5 + UP*1.5, radius=0.01).set_color(YELLOW))
                
        lights_dot.shuffle()

        self.play(ShowIncreasingSubsets(lights_dot))
        
        self.next_slide(loop=True)

        self.wait(3)

        self.next_slide(auto_next=True)

        wave_equation = Tex(r"$\hat{H}| \psi \rangle = E | \psi \rangle$")
        wave_equation.scale(0.8)

        interpretation = Tex(r"$|\psi|^2$ : probabilité de présence")
        interpretation.scale(0.8)

        text = Group(wave_equation, interpretation).arrange(buff=1, direction=DOWN)

        text.shift(RIGHT*3.5+DOWN*1.5)

        self.play(FadeIn(wave_equation))

        self.next_slide(loop=True)

        self.wait(3)

        self.next_slide(auto_next=True)

        self.play(FadeIn(interpretation))

        self.next_slide(loop=True)

        self.wait(3)

        self.next_slide()

        for m in self.mobjects:
            m.clear_updaters() 

        self.play(*[FadeOut(mob)for mob in self.mobjects])
        


    def create_a_wave(self, x, y):
        tracker = ValueTracker(0)
        wave = Arc(angle=PI/2, start_angle=-PI/4, radius=tracker.get_value()).set_color(YELLOW).shift(x*RIGHT+y*UP)
        wave.add_updater(lambda m: m.become(Arc(angle=PI/4, start_angle=-PI/8, radius=tracker.get_value())).set_color(YELLOW).shift(x*RIGHT+y*UP))
        tracker.add_updater(lambda m, dt: m.increment_value(dt) if tracker.get_value() < 3 else m.set_value(0))

        self.add(wave)
        self.add(tracker)

    def create_wavefront(self, x, y):
        self.create_a_wave(x, y)
        self.wait(.5)
        self.create_a_wave(x, y)
        self.wait(.5)
        self.create_a_wave(x, y)
        self.wait(.5)
        self.create_a_wave(x, y)
        self.wait(.5)
        self.create_a_wave(x, y)
        self.wait(.5)
        self.create_a_wave(x, y)
        self.wait(.5)

    def sym_wavefront(self, x, y):
        self.create_a_wave(x, y)
        self.create_a_wave(x, -y)
        self.wait(.5)
        self.create_a_wave(x, y)
        self.create_a_wave(x, -y)
        self.wait(.5)    
        self.create_a_wave(x, y)
        self.create_a_wave(x, -y)
        self.wait(.5)
        self.create_a_wave(x, y)
        self.create_a_wave(x, -y)
        self.wait(.5) 
        self.create_a_wave(x, y)
        self.create_a_wave(x, -y)
        self.wait(.5)
        self.create_a_wave(x, y)
        self.create_a_wave(x, -y)
        self.wait(.5)


class RevuePhysiqueSlide(Slide):
    def construct(self):
        title = Text("Notions de physique")
        title.shift(UP*3)
        title.scale(0.8)

        self.play(FadeIn(title))

        self.next_slide()

        schrodinger_equation = Tex(r"$\text{Équation de Schrödinger : } i \hbar \frac{\partial}{\partial t} |\psi\rangle = \hat{H} |\psi\rangle$")
        schrodinger_equation.shift(UP*2)

        self.play(Write(schrodinger_equation))

        self.next_slide()

        cadre_wave = RoundedRectangle(height=2.5, width=3.5, corner_radius=0.2, color=YELLOW)
        cadre_particle = RoundedRectangle(height=2.5, width=3.5, corner_radius=0.2, color=YELLOW)
        cadre_gpe = Group(cadre_wave, cadre_particle).arrange(buff=1)

        self.play(FadeIn(cadre_gpe))

        self.next_slide(auto_next=True)

        axes_sin = Axes(x_range=[-10, 10, 1], y_range=[-1, 1, 1], x_length=3, y_length=2).set_opacity(0)
        axes_particle = Axes(x_range=[-10, 10, 1], y_range=[-1, 1, 1], x_length=3, y_length=2)
        
        axes_gpe = Group(axes_sin, axes_particle).arrange(buff=1.5)

        tracker = ValueTracker(0)
        sin_graph = axes_sin.plot(lambda x: np.sin(x + tracker.get_value()), color=BLUE)
        sin_graph.add_updater(lambda m: m.become(axes_sin.plot(lambda x: np.sin(x + tracker.get_value()), color=BLUE)))

        self.play(Write(sin_graph))

        self.next_slide(loop=True)

        self.add(tracker)
        tracker.add_updater(lambda mobject, dt: mobject.increment_value(dt))
        self.wait(2*PI)

        self.next_slide(auto_next=True)

        particule = Dot(axes_particle.c2p(-5, 0), color=GREEN)
        momentum = Arrow(particule.get_center(), particule.get_center() + RIGHT*2, color=RED)
        momentum_label = Tex(r"$$\vec{p}$$")
        momentum_label.next_to(momentum, DOWN)

        particule_gpe = Group(particule, momentum, momentum_label)

        self.play(FadeIn(particule_gpe))

        self.next_slide(loop=True)

        self.wait(2*PI)

        self.next_slide(auto_next=True)

        sin_graph.set_z_index(1)
        axes_sin.set_z_index(-1)

        all_anim_phy = Group(axes_sin, cadre_gpe, particule_gpe)

        all_anim_phy.generate_target()
        all_anim_phy.target.shift(RIGHT*2)
        self.play(MoveToTarget(all_anim_phy))

        einstein_debroglie = Tex(r"$$E = h f$$ $$p = \frac{h}{\lambda}$$")
        einstein_debroglie.shift(LEFT*4)

        self.play(Write(einstein_debroglie))

        self.next_slide(loop=True)

        self.wait(2*PI)

        self.next_slide(auto_next=True)

        born_interpretation = Tex(r"$$\text{Interprétation de Born : } \int_{-\infty}^{+\infty} |\psi(x)|^2 dx = 1$$")
        born_interpretation.shift(DOWN*2)

        self.play(FadeIn(born_interpretation))

        self.next_slide(loop=True)

        self.wait(2*PI)

        self.next_slide()

        tracker.clear_updaters()
        sin_graph.clear_updaters()
        self.play(FadeOut(title, schrodinger_equation, einstein_debroglie, born_interpretation, all_anim_phy, sin_graph))


class LanguageMathSlide(Slide):
    def construct(self):
        title = Text("Language mathématique")
        title.shift(UP*3)

        self.play(FadeIn(title))

        self.next_slide()

        subt_info = Text("Informatique :")
        order_l1 = Tex(r"Notation : ", r"$\mathcal{O}(f(n))$", r" $\rightarrow$ ordre de grandeur")
        order_l1.set_color_by_tex("ordre de grandeur", GREEN)
        order_l1.set_color_by_tex("f(n)", YELLOW)
        order_l2 = Tex(r"example :", r" $\mathcal{O}(a \cdot n^3)$ ", r"$\sim$", r" $\mathcal{O}(n^3)$ ")
        order_l2.set_color_by_tex("n^3", BLUE)
        order_l2.set_color_by_tex("(a", RED)
        order_l3 = Tex(r"$\Rightarrow$ comparaison pour de grandes valeurs de $n$")

        info_gpe = Group(subt_info, order_l1, order_l2, order_l3).arrange(buff=0.5, direction=DOWN)

        self.play(FadeIn(subt_info))

        self.next_slide()

        self.play(Write(order_l1))

        self.next_slide()

        self.play(Write(order_l2))

        self.next_slide()

        self.play(Write(order_l3))

        self.next_slide()

        self.play(Transform(info_gpe, (info_gpe.scale(0.5)).shift(UP*1.5)))

        dem_line = Line(LEFT*5, RIGHT*5)
        self.play(Write(dem_line.shift(UP*0.4)))

        self.next_slide()

        subt_phy = Text('Notation de Dirac, dite "Bra-Ket" :')

        ket = Tex(r"$$|\psi\rangle = \begin{pmatrix} \psi_1 \\ \psi_2 \end{pmatrix}$$")
        matrix = Tex(r"$$M = \begin{pmatrix} m_1 & m_2 \\ m_3 & m_4 \end{pmatrix}$$")
        tensor = Tex(r"$$|a\rangle \otimes |b\rangle = \begin{pmatrix} a_1 b_1 \\ a_1 b_2 \\ a_2 b_1 \\ a_2 b_2 \end{pmatrix}$$")
        matrix_ket = Tex(r"$$M |\psi\rangle = \begin{pmatrix} m_1 \psi_1 + m_2 \psi_2 \\ m_3 \psi_1 + m_4 \psi_2 \end{pmatrix}$$")
        dot_expl = Tex(r"$$ \langle \psi | = (\psi_1^* \quad \psi_2^*) \rightarrow \langle \psi | \psi \rangle = |\psi|^2$$")

        gpe_m_et_ket = Group(ket, matrix, tensor, matrix_ket).arrange_in_grid(buff=0.5)
        gpe_phy = Group(subt_phy, gpe_m_et_ket, dot_expl).arrange(buff=0.5, direction=DOWN)

        gpe_phy.scale(0.5)
        gpe_phy.shift(DOWN*1.5)

        self.play(FadeIn(subt_phy))

        self.next_slide()

        self.play(Write(ket))

        self.next_slide()

        self.play(Write(matrix))

        self.next_slide()

        self.play(Write(tensor))

        self.next_slide()

        self.play(Write(matrix_ket))

        self.next_slide()

        self.play(Write(dot_expl))

        self.next_slide()

        self.play(FadeOut(title, info_gpe, dem_line, gpe_phy))


class QubitSlide(ThreeDSlide):
    def construct(self):
        title = Tex("Des qubits")
        title.move_to([0, 3, 0])

        self.add_fixed_in_frame_mobjects(title)
        self.play(FadeIn(title))

        self.next_slide()

        axes = ThreeDAxes(
            x_range = (-1, 1, 1), 
            y_range = (-1, 1, 1), 
            z_range = (-1, 1, 1), 
            x_length = 7, 
            y_length = 7, 
            z_length = 7
        )

        self.move_camera(
            phi = 70*DEGREES, 
            theta = -35*DEGREES,
            distance = 1
        )

        arrow = Arrow3D(
            start = np.array([0, 0, 0]),
            end = np.array([1, 1, 1]),
            resolution = 8, 
            color = BLUE, 
        )

        
        bloch_zero = Arrow3D(
            start = np.array([0, 0, 0]),
            end = np.array([0, 0, 2]),
            resolution = 8, 
            color = BLUE, 
        )

        bloch_one = Arrow3D(
            start = np.array([0, 0, 0]),
            end = np.array([0, 0, -2]),
            resolution = 8, 
            color = BLUE, 
        )

        bloch_plus = Arrow3D(
            start = np.array([0, 0, 0]),
            end = np.array([2, 0, 0]),
            resolution = 8, 
            color = BLUE, 
        )

        phi, theta, focal_distance, gamma, distance_to_origin = self.camera.get_value_trackers()

        sphere = Surface(
            lambda u, v: np.array([
            2 * np.cos(u) * np.cos(v),
            2 * np.cos(u) * np.sin(v),
            2 * np.sin(u)
            ]), v_range=[0, TAU], u_range=[-PI / 2, PI / 2],
            checkerboard_colors=[RED, PURPLE], resolution=(16, 32), fill_opacity = 0.1
        )
        sphere.set_color(YELLOW)

        one = MathTex(r"|1\rangle").move_to([0.6, 0, -3]).scale(0.6)
        zero = MathTex(r"|0\rangle").move_to([0.6, 0, 3]).scale(0.6)

        self.add_fixed_orientation_mobjects(one)
        self.add_fixed_orientation_mobjects(zero)

        all_bloch = Group(axes, arrow, sphere, one, zero, bloch_zero, bloch_one, bloch_plus)
        all_bloch.scale(0.6)
        all_bloch.move_to([-5, -2, 0])

        self.play(Create(axes), Write(one), Write(zero), run_time = 2)

        self.play(GrowFromPoint(arrow, [-5, -2, 0]), Create(sphere), run_time = 4)

        qubit_on_bloch = Tex(r"$\left| \psi \right\rangle = \alpha \left| 0 \right\rangle + \beta \left| 1 \right\rangle, (\alpha, \beta) \in \mathbb{C}^2$")
        qubit_on_bloch.set_color(BLUE)
        qubit_on_bloch.scale(0.6)
        qubit_on_bloch.move_to([-4, -1, 0])

        self.add_fixed_in_frame_mobjects(qubit_on_bloch)

        self.play(Write(qubit_on_bloch))

        self.next_slide()

        qcircuit_temp = TexTemplate()
        qcircuit_temp.add_to_preamble(r"\usepackage[braket, qm]{qcircuit} \usepackage{graphicx}")

        circuit = Tex(r"""
        \scalebox{1.0}{
        \Qcircuit @C=1.0em @R=0.2em @!R { \\
        	 	\nghost{ } & \lstick{ } & \gate{\mathrm{X}} & \targ & \control\qw \barrier[0em]{1} & \qw & \meter & \qw & \qw & \qw\\
        	 	\nghost{ } & \lstick{ } & \gate{\mathrm{H}} & \ctrl{-1} & \ctrl{-1} & \qw & \qw & \meter & \qw & \qw\\
        	 	\nghost{ } & \lstick{ } & \lstick{/_{_{2}}} \cw & \cw & \cw & \cw & \dstick{_{_{\hspace{0.0em}0} } } \cw \ar @{<=} [-2,0] & \dstick{_{_{\hspace{0.0em}1} } } \cw \ar @{<=} [-1,0] & \cw & \cw\\
        \\ }}
        """, tex_template=qcircuit_temp)
        circuit.scale(0.8)
        circuit.move_to([3, 1, 0])

        self.add_fixed_in_frame_mobjects(circuit)
        self.play(FadeIn(circuit))

        self.next_slide()

        self.play(ReplacementTransform(arrow, bloch_one), run_time = 2)

        qubit_one = Tex(r"$\left| \psi \right\rangle = \left| 1 \right\rangle$")
        qubit_one.scale(0.6)
        qubit_one.move_to([-4, -2, 0])

        self.add_fixed_in_frame_mobjects(qubit_one)

        self.play(FadeIn(qubit_one))

        self.next_slide()

        box_x = Square(side_length=0.6, color=ORANGE).move_to([0.87, 1.6, 0])
        self.add_fixed_in_frame_mobjects(box_x)
        self.play(Write(box_x))

        self.play(ReplacementTransform(bloch_one, bloch_zero), run_time = 2)

        qubit_zero = Tex(r"$\mathrm{X} \left| 1 \right\rangle = \left| 0 \right\rangle$")
        qubit_zero.set_color(ORANGE)
        qubit_zero.scale(0.6)
        qubit_zero.move_to([-4, -3, 0])

        self.add_fixed_in_frame_mobjects(qubit_zero)

        self.play(FadeIn(qubit_zero))

        self.next_slide()

        box_h = Square(side_length=0.6, color=RED).move_to([0.87, 0.97, 0])
        self.add_fixed_in_frame_mobjects(box_h)
        self.play(Write(box_h))

        self.play(ReplacementTransform(bloch_zero, bloch_plus), run_time = 2)

        qubit_plus = Tex(r"$\mathrm{H} \left| 0 \right\rangle = \left| + \right\rangle = \frac{1}{\sqrt{2}} \left| 0 \right\rangle + \frac{1}{\sqrt{2}} \left| 1 \right\rangle$")
        qubit_plus.set_color(RED)
        qubit_plus.scale(0.6)
        qubit_plus.move_to([4, -1, 0])

        self.add_fixed_in_frame_mobjects(qubit_plus)

        self.play(FadeIn(qubit_plus))

        self.next_slide()

        box_cx = Rectangle(height=1.2, width=0.6, color=GREEN).move_to([1.7, 1.3, 0])
        self.add_fixed_in_frame_mobjects(box_cx)
        self.play(Write(box_cx))

        qubit_int = Tex(r"$\mathrm{CX} (\left| + \right\rangle_c \otimes \left| 0 \right\rangle_t)  = \frac{1}{\sqrt{2}} \left| 00 \right\rangle + \frac{1}{\sqrt{2}} \left| 11 \right\rangle$")
        qubit_int.set_color(GREEN)
        qubit_int.scale(0.6)
        qubit_int.move_to([4, -2, 0])

        self.add_fixed_in_frame_mobjects(qubit_int)

        self.play(FadeIn(qubit_int))

        self.next_slide()

        box_cz = Rectangle(height=1.2, width=0.6, color=PURPLE).move_to([2.4, 1.3, 0])
        self.add_fixed_in_frame_mobjects(box_cz)
        self.play(Write(box_cz))

        qubit_cz = Tex(r"$\mathrm{CZ} (\frac{1}{\sqrt{2}} \left| 00 \right\rangle + \frac{1}{\sqrt{2}} \left| 11 \right\rangle) = \frac{1}{\sqrt{2}} \left| 00 \right\rangle - \frac{1}{\sqrt{2}} \left| 11 \right\rangle$")
        qubit_cz.set_color(PURPLE)
        qubit_cz.scale(0.6)
        qubit_cz.move_to([4, -3, 0])

        self.add_fixed_in_frame_mobjects(qubit_cz)

        self.play(FadeIn(qubit_cz))

        self.next_slide()

        box_mes = Rectangle(height=1.4, width=2, color=YELLOW).move_to([4.1, 1.3, 0])
        self.add_fixed_in_frame_mobjects(box_mes)
        self.play(Write(box_mes))

        mesure = Tex(r"Mesure : 00 ou 11")
        mesure.set_color(YELLOW)
        mesure.scale(0.6)
        mesure.move_to([4.1, 2.4, 0])

        self.add_fixed_in_frame_mobjects(mesure)

        self.play(FadeIn(mesure))

        self.next_slide()

        self.play(*[FadeOut(mob)for mob in self.mobjects])        


class ExemplesAlgoSlide(Slide):
    def construct(self):
        title = Text("Exemples d'algorithmes")
        title.shift(UP*3)

        self.play(FadeIn(title))

        self.next_slide()

        qcircuit_temp = TexTemplate()
        qcircuit_temp.add_to_preamble(r"\usepackage[braket, qm]{qcircuit} \usepackage{graphicx}")

        dj_circuit = Tex(r"""
        \scalebox{1.0}{
        \Qcircuit @C=1.0em @R=0.2em @!R { \\
        	 	\nghost{ {input} :  } & \lstick{ {input} :  } & \qw \barrier[0em]{1} & \qw & \gate{\mathrm{H} } & \multigate{1}{\mathrm{Oracle} }_<<<{0} & \gate{\mathrm{H} } & \meter & \qw & \qw\\
        	 	\nghost{ {output} :  } & \lstick{ {output} :  } & \gate{\mathrm{X} } & \qw & \gate{\mathrm{H} } & \ghost{\mathrm{Oracle} }_<<<{1} & \qw & \qw & \qw & \qw\\
        	 	\nghost{\mathrm{ {measure} :  } } & \lstick{\mathrm{ {measure} :  } } & \lstick{/_{_{1}} } \cw & \cw & \cw & \cw & \cw & \dstick{_{_{\hspace{0.0em}0} } } \cw \ar @{<=} [-2,0] & \cw & \cw\\
        \\ }}
        """, tex_template=qcircuit_temp)
        
        shor_circuit = Tex(r"""
        \scalebox{1.0}{
        \Qcircuit @C=1.0em @R=0.2em @!R { \\
        	 	\nghost{ {q0}_{0} :  } & \lstick{ {q0}_{0} :  } & \gate{\mathrm{H}} & \ctrl{8} & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \multigate{7}{\mathrm{QFT†}}_<<<{0} & \meter & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {q0}_{1} :  } & \lstick{ {q0}_{1} :  } & \gate{\mathrm{H}} & \qw & \ctrl{7} & \qw & \qw & \qw & \qw & \qw & \qw & \ghost{\mathrm{QFT†}}_<<<{1} & \qw & \meter & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {q0}_{2} :  } & \lstick{ {q0}_{2} :  } & \gate{\mathrm{H}} & \qw & \qw & \ctrl{6} & \qw & \qw & \qw & \qw & \qw & \ghost{\mathrm{QFT†}}_<<<{2} & \qw & \qw & \meter & \qw & \qw & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {q0}_{3} :  } & \lstick{ {q0}_{3} :  } & \gate{\mathrm{H}} & \qw & \qw & \qw & \ctrl{5} & \qw & \qw & \qw & \qw & \ghost{\mathrm{QFT†}}_<<<{3} & \qw & \qw & \qw & \meter & \qw & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {q0}_{4} :  } & \lstick{ {q0}_{4} :  } & \gate{\mathrm{H}} & \qw & \qw & \qw & \qw & \ctrl{4} & \qw & \qw & \qw & \ghost{\mathrm{QFT†}}_<<<{4} & \qw & \qw & \qw & \qw & \meter & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {q0}_{5} :  } & \lstick{ {q0}_{5} :  } & \gate{\mathrm{H}} & \qw & \qw & \qw & \qw & \qw & \ctrl{3} & \qw & \qw & \ghost{\mathrm{QFT†}}_<<<{5} & \qw & \qw & \qw & \qw & \qw & \meter & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {q0}_{6} :  } & \lstick{ {q0}_{6} :  } & \gate{\mathrm{H}} & \qw & \qw & \qw & \qw & \qw & \qw & \ctrl{2} & \qw & \ghost{\mathrm{QFT†}}_<<<{6} & \qw & \qw & \qw & \qw & \qw & \qw & \meter & \qw & \qw & \qw\\
        	 	\nghost{ {q0}_{7} :  } & \lstick{ {q0}_{7} :  } & \gate{\mathrm{H}} & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \ctrl{1} & \ghost{\mathrm{QFT†}}_<<<{7} & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \meter & \qw & \qw\\
        	 	\nghost{ {q0}_{8} :  } & \lstick{ {q0}_{8} :  } & \gate{\mathrm{X}} & \multigate{3}{\mathrm{7\string^1\,mod\,15}}_<<<{0} & \multigate{3}{\mathrm{7\string^2\,mod\,15}}_<<<{0} & \multigate{3}{\mathrm{7\string^4\,mod\,15}}_<<<{0} & \multigate{3}{\mathrm{7\string^8\,mod\,15}}_<<<{0} & \multigate{3}{\mathrm{7\string^16\,mod\,15}}_<<<{0} & \multigate{3}{\mathrm{7\string^32\,mod\,15}}_<<<{0} & \multigate{3}{\mathrm{7\string^64\,mod\,15}}_<<<{0} & \multigate{3}{\mathrm{7\string^128\,mod\,15}}_<<<{0} & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {q0}_{9} :  } & \lstick{ {q0}_{9} :  } & \qw & \ghost{\mathrm{7\string^1\,mod\,15}}_<<<{1} & \ghost{\mathrm{7\string^2\,mod\,15}}_<<<{1} & \ghost{\mathrm{7\string^4\,mod\,15}}_<<<{1} & \ghost{\mathrm{7\string^8\,mod\,15}}_<<<{1} & \ghost{\mathrm{7\string^16\,mod\,15}}_<<<{1} & \ghost{\mathrm{7\string^32\,mod\,15}}_<<<{1} & \ghost{\mathrm{7\string^64\,mod\,15}}_<<<{1} & \ghost{\mathrm{7\string^128\,mod\,15}}_<<<{1} & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {q0}_{10} :  } & \lstick{ {q0}_{10} :  } & \qw & \ghost{\mathrm{7\string^1\,mod\,15}}_<<<{2} & \ghost{\mathrm{7\string^2\,mod\,15}}_<<<{2} & \ghost{\mathrm{7\string^4\,mod\,15}}_<<<{2} & \ghost{\mathrm{7\string^8\,mod\,15}}_<<<{2} & \ghost{\mathrm{7\string^16\,mod\,15}}_<<<{2} & \ghost{\mathrm{7\string^32\,mod\,15}}_<<<{2} & \ghost{\mathrm{7\string^64\,mod\,15}}_<<<{2} & \ghost{\mathrm{7\string^128\,mod\,15}}_<<<{2} & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {q0}_{11} :  } & \lstick{ {q0}_{11} :  } & \qw & \ghost{\mathrm{7\string^1\,mod\,15}}_<<<{3} & \ghost{\mathrm{7\string^2\,mod\,15}}_<<<{3} & \ghost{\mathrm{7\string^4\,mod\,15}}_<<<{3} & \ghost{\mathrm{7\string^8\,mod\,15}}_<<<{3} & \ghost{\mathrm{7\string^16\,mod\,15}}_<<<{3} & \ghost{\mathrm{7\string^32\,mod\,15}}_<<<{3} & \ghost{\mathrm{7\string^64\,mod\,15}}_<<<{3} & \ghost{\mathrm{7\string^128\,mod\,15}}_<<<{3} & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{\mathrm{ {c0} :  } } & \lstick{\mathrm{ {c0} :  } } & \lstick{/_{_{8}} } \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \dstick{_{_{\hspace{0.0em}0} } } \cw \ar @{<=} [-12,0] & \dstick{_{_{\hspace{0.0em}1} } } \cw \ar @{<=} [-11,0] & \dstick{_{_{\hspace{0.0em}2} } } \cw \ar @{<=} [-10,0] & \dstick{_{_{\hspace{0.0em}3} } } \cw \ar @{<=} [-9,0] & \dstick{_{_{\hspace{0.0em}4} } } \cw \ar @{<=} [-8,0] & \dstick{_{_{\hspace{0.0em}5} } } \cw \ar @{<=} [-7,0] & \dstick{_{_{\hspace{0.0em}6} } } \cw \ar @{<=} [-6,0] & \dstick{_{_{\hspace{0.0em}7} } } \cw \ar @{<=} [-5,0] & \cw & \cw\\
        \\ }}
        """, tex_template=qcircuit_temp)

        key_circuit = Tex(r"""
        \scalebox{1.0}{
        \Qcircuit @C=1.0em @R=0.2em @!R { \\
        	 	\nghost{ {q} :  } & \lstick{ {q} :  } & \gate{\mathrm{H} } \barrier[0em]{0} & \qw & \gate{\mathrm{H} } & \meter & \qw & \qw\\
        	 	\nghost{\mathrm{ {c} :  } } & \lstick{\mathrm{ {c} :  } } & \lstick{/_{_{1}}} \cw & \cw & \cw & \dstick{_{_{\hspace{0.0em}0} } } \cw \ar @{<=} [-1,0] & \cw & \cw\\
        \\ } }

        \scalebox{1.0}{
        \Qcircuit @C=1.0em @R=0.2em @!R { \\
        	 	\nghost{ {q} :  } & \lstick{ {q} :  } & \gate{\mathrm{H} } \barrier[0em]{0} & \qw & \meter \barrier[0em]{0} & \qw & \gate{\mathrm{H} } & \meter & \qw & \qw\\
        	 	\nghost{\mathrm{ {c} :  } } & \lstick{\mathrm{ {c} :  } } & \lstick{/_{_{1}} } \cw & \cw & \dstick{_{_{\hspace{0.0em}0} } } \cw \ar @{<=} [-1,0] & \cw & \cw & \dstick{_{_{\hspace{0.0em}0} } } \cw \ar @{<=} [-1,0] & \cw & \cw\\
        \\ } }
        """, tex_template=qcircuit_temp)

        grov_circuit = Tex(r"""
        \scalebox{1.0}{
        \Qcircuit @C=1.0em @R=0.2em @!R { \\
        	 	\nghost{ {q}_{0} :  } & \lstick{ {q}_{0} :  } & \gate{\mathrm{H} } \barrier[0em]{3} & \qw & \ctrl{1} & \ctrl{2} & \ctrl{3} \barrier[0em]{3} & \qw & \gate{\mathrm{H} } & \gate{\mathrm{X} } & \ctrl{1} & \gate{\mathrm{X} } & \gate{\mathrm{H} } \barrier[0em]{3} & \qw \barrier[0em]{3} & \qw & \meter & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {q}_{1} :  } & \lstick{ {q}_{1} :  } & \gate{\mathrm{H} } & \qw & \control\qw & \qw & \qw & \qw & \gate{\mathrm{H} } & \gate{\mathrm{X} } & \ctrl{1} & \gate{\mathrm{X} } & \gate{\mathrm{H} } & \qw & \qw & \qw & \meter & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {q}_{2} :  } & \lstick{ {q}_{2} :  } & \gate{\mathrm{H} } & \qw & \qw & \control\qw & \qw & \qw & \gate{\mathrm{H} } & \gate{\mathrm{X} } & \ctrl{1} & \gate{\mathrm{X} } & \gate{\mathrm{H} } & \qw & \qw & \qw & \qw & \meter & \qw & \qw & \qw\\
        	 	\nghost{ {q}_{3} :  } & \lstick{ {q}_{3} :  } & \gate{\mathrm{H} } & \qw & \qw & \qw & \control\qw & \qw & \gate{\mathrm{H} } & \gate{\mathrm{X} } & \control\qw & \gate{\mathrm{X} } & \gate{\mathrm{H} } & \qw & \qw & \qw & \qw & \qw & \meter & \qw & \qw\\
        	 	\nghost{\mathrm{ {meas} :  } } & \lstick{\mathrm{ {meas} :  } } & \lstick{/_{_{4} } } \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \dstick{_{_{\hspace{0.0em}0} } } \cw \ar @{<=} [-4,0] & \dstick{_{_{\hspace{0.0em}1} } } \cw \ar @{<=} [-3,0] & \dstick{_{_{\hspace{0.0em}2} } } \cw \ar @{<=} [-2,0] & \dstick{_{_{\hspace{0.0em}3} } } \cw \ar @{<=} [-1,0] & \cw & \cw\\
        \\ } }
        """, tex_template=qcircuit_temp)

        iter_circuit = Tex(r"""
        \scalebox{1.0}{
        \Qcircuit @C=1.0em @R=0.2em @!R { \\
        	 	\nghost{ {q}_{0} :  } & \lstick{ {q}_{0} :  } & \gate{\mathrm{H}} & \control \qw & \dstick{\hspace{2.0em}\mathrm{P}\,(\mathrm{\frac{\pi}{4}})} \qw & \qw & \qw & \control \qw & \dstick{\hspace{2.0em}\mathrm{P}\,(\mathrm{\frac{\pi}{4}})} \qw & \qw & \qw & \control \qw & \dstick{\hspace{2.0em}\mathrm{P}\,(\mathrm{\frac{\pi}{4}})} \qw & \qw & \qw & \control \qw & \dstick{\hspace{2.0em}\mathrm{P}\,(\mathrm{\frac{\pi}{4}})} \qw & \qw & \qw & \gate{\mathrm{H}} & \meter \barrier[0em]{1} & \qw & \gate{\mathrm{\left|0\right\rangle}} & \gate{\mathrm{H}} & \control \qw & \dstick{\hspace{2.0em}\mathrm{P}\,(\mathrm{\frac{\pi}{4}})} \qw & \qw & \qw & \control \qw & \dstick{\hspace{2.0em}\mathrm{P}\,(\mathrm{\frac{\pi}{4}})} \qw & \qw & \qw & \gate{\mathrm{P}\,(\mathrm{\frac{-\pi}{2}})} & \gate{\mathrm{H}} & \meter \barrier[0em]{1} & \qw & \gate{\mathrm{\left|0\right\rangle}} & \gate{\mathrm{H}} & \control \qw & \dstick{\hspace{2.0em}\mathrm{P}\,(\mathrm{\frac{\pi}{4}})} \qw & \qw & \qw & \gate{\mathrm{P}\,(\mathrm{\frac{-\pi}{2}})} & \gate{\mathrm{P}\,(\mathrm{\frac{-\pi}{4}})} & \gate{\mathrm{H}} & \meter & \qw & \qw\\
        	 	\nghost{ {q}_{1} :  } & \lstick{ {q}_{1} :  } & \gate{\mathrm{X}} & \ctrl{-1} & \qw & \qw & \qw & \ctrl{-1} & \qw & \qw & \qw & \ctrl{-1} & \qw & \qw & \qw & \ctrl{-1} & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \ctrl{-1} & \qw & \qw & \qw & \ctrl{-1} & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \ctrl{-1} & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{\mathrm{ {c4} :  } } & \lstick{\mathrm{ {c4} :  } } & \lstick{/_{_{3} } } \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \dstick{_{_{\hspace{0.0em}0} } } \cw \ar @{<=} [-2,0] & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \control \cw^(0.0){^{\mathtt{c4_0=0x1} } } \cwx[-2] & \cw & \dstick{_{_{\hspace{0.0em}1} } } \cw \ar @{<=} [-2,0] & \cw & \cw & \cw & \cw & \cw & \cw & \cw & \control \cw^(0.0){^{\mathtt{c4_1=0x1} } } \cwx[-2] & \control \cw^(0.0){^{\mathtt{c4_0=0x1} } } \cwx[-2] & \cw & \dstick{_{_{\hspace{0.0em}2} } } \cw \ar @{<=} [-2,0] & \cw & \cw\\
        \\ }}
        """, tex_template=qcircuit_temp)

        circuit_with_title = []
        for i in [[Text("Deutsch-Josza"), dj_circuit.scale(2)], [Text("Shor"), shor_circuit], [Text("QKD"), key_circuit.scale(2)], [Text("Grover"), grov_circuit.scale(1.5)], [Text("IPE"), iter_circuit.scale(1.5)]]:
            circuit_with_title.append(Group(i[0].scale(0.3), i[1].scale(0.1)).arrange(buff=0.5, direction=DOWN))

        up_part = Group(*circuit_with_title[0:3]).arrange(buff=1, direction=RIGHT)
        down_part = Group(*circuit_with_title[3:5]).arrange(buff=1, direction=RIGHT)

        all_circ = Group(up_part, down_part).arrange(buff=1, direction=DOWN)

        color_list = [YELLOW, BLUE, RED, GREEN, PURPLE]
        for i in range(len(circuit_with_title)):
            circuit_with_title[i].set_color(color_list[i])

        self.play(Write(circuit_with_title[0][0]))
        self.play(Write(circuit_with_title[0][1]))

        self.next_slide()

        self.play(Write(circuit_with_title[1][0]))
        self.play(Write(circuit_with_title[1][1]))

        self.next_slide()

        self.play(Write(circuit_with_title[2][0]))
        self.play(Write(circuit_with_title[2][1]))

        self.next_slide()

        self.play(Write(circuit_with_title[3][0]))
        self.play(Write(circuit_with_title[3][1]))

        self.next_slide()

        self.play(Write(circuit_with_title[4][0]))
        self.play(Write(circuit_with_title[4][1]))

        self.next_slide()

        self.play(FadeOut(title, *circuit_with_title))


class DeutschSlide(Slide):
    def construct(self):
        title = Text("Algorithme de Deutsch")
        title.shift(UP*3)

        self.play(FadeIn(title))

        self.next_slide()

        qcircuit_temp = TexTemplate()
        qcircuit_temp.add_to_preamble(r"\usepackage[braket, qm]{qcircuit} \usepackage{graphicx}")

        dj_circuit = Tex(r"""
        \scalebox{1.0}{
        \Qcircuit @C=1.0em @R=0.2em @!R { \\
        	 	\nghost{ {input} :  } & \lstick{ {input} :  } & \qw \barrier[0em]{1} & \qw & \gate{\mathrm{H} } & \multigate{1}{\mathrm{Oracle} }_<<<{0} & \gate{\mathrm{H} } & \meter & \qw & \qw\\
        	 	\nghost{ {output} :  } & \lstick{ {output} :  } & \gate{\mathrm{X} } & \qw & \gate{\mathrm{H} } & \ghost{\mathrm{Oracle} }_<<<{1} & \qw & \qw & \qw & \qw\\
        	 	\nghost{\mathrm{ {measure} :  } } & \lstick{\mathrm{ {measure} :  } } & \lstick{/_{_{1}} } \cw & \cw & \cw & \cw & \cw & \dstick{_{_{\hspace{0.0em}0} } } \cw \ar @{<=} [-2,0] & \cw & \cw\\
        \\ }}
        """, tex_template=qcircuit_temp)
        dj_circuit.scale(0.5)
        dj_circuit.shift(UP*0.8)

        problem = Tex(r"Problème : Soit $f : \{0, 1\} \rightarrow \{0, 1\}$, tel que $f$ est soit constante ($f(0)=f(1)$) soit équilibrée ($f(0)\neq f(1)$). Trouver si $f$ est constante ou équilibrée.", color=YELLOW)
        problem.scale(0.5)
        problem.shift(UP*2)

        self.play(FadeIn(problem))

        self.next_slide()

        self.play(Write(dj_circuit))

        self.next_slide()

        init = Tex(r"Initialisation : $\left| 0 0 \right\rangle \rightarrow \left| 0 1 \right\rangle$")
        appl_h = Tex(r"Application de $\mathrm{H}^{\otimes 2}$ : $$\left| 0 1 \right\rangle \rightarrow \frac{1}{\sqrt{2}} (\left| 0 \right\rangle + \left| 1 \right\rangle) \cdot \frac{1}{\sqrt{2}} (\left| 0 \right\rangle - \left| 1 \right\rangle) = \frac{1}{2} (\left| 0 \right\rangle (\left| 0 \right\rangle - \left| 1 \right\rangle) + \left| 1 \right\rangle (\left| 0 \right\rangle - \left| 1 \right\rangle)) $$")
        oracle = Tex(r"Oracle : $\left| x \right\rangle \left| y \right\rangle = \left| x \right\rangle \left| f(x) \oplus y \right\rangle$")
        appl_oracle = Tex(r"Application de l'oracle : $$\frac{1}{2} (\left| 0 \right\rangle (\left| f(0) \oplus 0 \right\rangle - \left| f(0) \oplus 1 \right\rangle) + \left| 1 \right\rangle (\left| f(1) \oplus 0 \right\rangle - \left| f(1) \oplus 1 \right\rangle)) = \pm \frac{1}{2} (\left| 0 \right\rangle \pm \left| 1 \right\rangle) (\left| 0 \right\rangle - \left| 1 \right\rangle)$$")
        concl = Tex(r"Cas $\left| 0 \right\rangle + \left| 1 \right\rangle$ [$f$ constante] : après $\mathrm{H}$, cela passe à $\left| 0 \right\rangle$ ; cas $\left| 0 \right\rangle - \left| 1 \right\rangle$ [$f$ équilibrée] : après $\mathrm{H}$,  cela passe à $\left| 1 \right\rangle$", color=GREEN)

        expl = Group(init, appl_h, oracle, appl_oracle, concl).arrange(buff=0.5, direction=DOWN)
        expl.scale(0.4)
        expl.shift(DOWN*2)

        self.play(FadeIn(expl))

        self.next_slide()

        diff_class = Tex(r"Complexité classique : $\mathcal{O}(2^n)$, et quantique : $\mathcal{O}(1)$ !", color=RED)
        diff_class.scale(0.5)
        diff_class.shift(UP*0)

        self.play(FadeIn(diff_class))

        self.next_slide()

        self.play(*[FadeOut(mob)for mob in self.mobjects])


class GroverSlide(Slide):
    def construct(self):
        title = Text("Algorithme de Grover")
        title.shift(UP*3)

        self.play(FadeIn(title))

        self.next_slide()

        qcircuit_temp = TexTemplate()
        qcircuit_temp.add_to_preamble(r"\usepackage[braket, qm]{qcircuit} \usepackage{graphicx} \usepackage{physics}")

        grov_circuit = Tex(r"""
        \scalebox{1.0}{
            \Qcircuit @C=1.0em @R=0.2em @!R { \\
            \nghost{\ket{0}} & \lstick{\ket{0}} & {/^n} \qw & \gate{\mathrm{H^{\otimes n}}} & \gate{Oracle} & \gate{Diffuseur} & \meter & \qw \gategroup{2}{5}{2}{6}{.7em}{--} \\
            & & & &  \mbox{$\qquad \qquad \quad$ Répéter $\order{\sqrt{n}}$ fois} & & & \\
            \\ }}
        """, tex_template=qcircuit_temp)

        grov_circuit.scale(0.8)
        grov_circuit.shift(DOWN*0.8)

        doc = fitz.open(r"step1-geo.pdf")

        page = doc.load_page(0)
        svg = page.get_svg_image()
        
        file = open(r"media/step1.svg", "w")
        file.write(svg)
        file.close()

        step1 = SVGMobject(r"media/step1.svg").set_color(YELLOW).scale_to_fit_width(4)

        doc = fitz.open(r"step2-geo.pdf")

        page = doc.load_page(0)
        svg = page.get_svg_image()
        
        file = open(r"media/step2.svg", "w")
        file.write(svg)
        file.close()

        step2 = SVGMobject(r"media/step2.svg").set_color(GREEN).scale_to_fit_width(4)

        doc = fitz.open(r"step3-geo.pdf")

        page = doc.load_page(0)
        svg = page.get_svg_image()
        
        file = open(r"media/step3.svg", "w")
        file.write(svg)
        file.close()

        step3 = SVGMobject(r"media/step3.svg").set_color(BLUE).scale_to_fit_width(4)

        schem = Group(step1, step2, step3).arrange(buff=0.5, direction=RIGHT).shift(UP*1.5)

        self.play(Write(grov_circuit))

        self.next_slide()

        step1_cir = Rectangle(height=0.7, width=1.1, color=YELLOW).shift(LEFT*2.28, DOWN*0.6)

        self.play(Write(step1_cir))
        self.play(FadeIn(step1))

        self.next_slide()

        step2_cir = Rectangle(height=0.7, width=1.6, color=GREEN).shift(LEFT*0.65, DOWN*0.6)

        self.play(Write(step2_cir))
        self.play(FadeIn(step2))

        self.next_slide()

        step3_cir = Rectangle(height=0.7, width=2.35, color=BLUE).shift(RIGHT*1.55, DOWN*0.6)

        self.play(Write(step3_cir))
        self.play(FadeIn(step3))

        self.next_slide()

        avant1 = Tex(r"Classique : $\mathcal{O}(N)$, et quantique : $\mathcal{O}(\sqrt{N})$ !", color=RED).shift(DOWN*2)
        avant2 = Tex(r"Classique : $\mathcal{O}(2^n)$, et quantique : $\mathcal{O}(2^{n/2})$ !", color=RED).shift(DOWN*2)

        self.play(FadeIn(avant1))

        self.next_slide()

        self.play(Transform(avant1, avant2))

        self.next_slide()

        self.play(*[FadeOut(mob)for mob in self.mobjects])


class SudokuSlide(Slide):
    def construct(self):
        title = Text("Algorithme de Grover - Sudoku")
        title.shift(UP*3)

        self.play(FadeIn(title))

        self.next_slide()

        empty_sudoku = Table([["1", "", ""],
                              ["2", "", "1"],
                              [ "", "1", ""]], include_outer_lines=True)
        empty_sudoku.scale(0.5)
        empty_sudoku.shift(LEFT*4+UP*1.5)

        fill_sudoku = Table([["1", "2", "3"],
                             ["2", "3", "1"],
                             ["3", "1", "2"]], include_outer_lines=True)
        fill_sudoku.add_highlighted_cell((1,2), color=GREEN)
        fill_sudoku.add_highlighted_cell((1,3), color=GREEN)
        fill_sudoku.add_highlighted_cell((2,2), color=GREEN)
        fill_sudoku.add_highlighted_cell((3,1), color=GREEN)
        fill_sudoku.add_highlighted_cell((3,3), color=GREEN)
        fill_sudoku.scale(0.5)
        fill_sudoku.shift(LEFT*4+UP*1.5)

        qcircuit_temp = TexTemplate()
        qcircuit_temp.add_to_preamble(r"\usepackage[braket, qm]{qcircuit} \usepackage{graphicx} \usepackage{physics}")

        oracle = Tex(r"""
        \scalebox{1.0}{
        \Qcircuit @C=1.0em @R=0.2em @!R { \\
        	 	\nghost{ {input}_{0} :  } & \lstick{ {input}_{0} :  } & \ctrl{5} & \qw & \qw & \ctrl{6} & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \ctrl{5} & \qw & \qw & \ctrl{6} & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {input}_{1} :  } & \lstick{ {input}_{1} :  } & \qw & \qw & \ctrl{4} & \qw & \qw & \qw & \ctrl{6} & \qw & \qw & \qw & \qw & \qw & \qw & \ctrl{4} & \qw & \qw & \qw & \ctrl{6} & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {input}_{2} :  } & \lstick{ {input}_{2} :  } & \qw & \ctrl{6} & \qw & \qw & \ctrl{7} & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \ctrl{6} & \qw & \qw & \ctrl{7} & \qw & \qw & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {input}_{3} :  } & \lstick{ {input}_{3} :  } & \qw & \qw & \qw & \qw & \qw & \ctrl{3} & \qw & \qw & \ctrl{7} & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \ctrl{3} & \qw & \qw & \ctrl{7} & \qw & \qw & \qw\\
        	 	\nghost{ {input}_{4} :  } & \lstick{ {input}_{4} :  } & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \ctrl{5} & \qw & \ctrl{6} & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \ctrl{5} & \qw & \ctrl{6} & \qw & \qw\\
        	 	\nghost{ {cond}_{0} :  } & \lstick{ {cond}_{0} :  } & \targ & \qw & \targ & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \ctrl{1} & \targ & \qw & \targ & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {cond}_{1} :  } & \lstick{ {cond}_{1} :  } & \qw & \qw & \qw & \targ & \qw & \targ & \qw & \qw & \qw & \qw & \ctrl{1} & \qw & \qw & \qw & \targ & \qw & \targ & \qw & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {cond}_{2} :  } & \lstick{ {cond}_{2} :  } & \qw & \qw & \qw & \qw & \qw & \qw & \targ & \qw & \qw & \qw & \ctrl{1} & \qw & \qw & \qw & \qw & \qw & \qw & \targ & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {cond}_{3} :  } & \lstick{ {cond}_{3} :  } & \qw & \targ & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \ctrl{1} & \qw & \targ & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {cond}_{4} :  } & \lstick{ {cond}_{4} :  } & \qw & \qw & \qw & \qw & \targ & \qw & \qw & \targ & \qw & \qw & \ctrl{1} & \qw & \qw & \qw & \qw & \targ & \qw & \qw & \targ & \qw & \qw & \qw & \qw\\
        	 	\nghost{ {cond}_{5} :  } & \lstick{ {cond}_{5} :  } & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \targ & \targ & \ctrl{1} & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \targ & \targ & \qw & \qw\\
        	 	\nghost{ {kick} :  } & \lstick{ {kick} :  } & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \targ & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw & \qw\\
        \\ } }
        """, tex_template=qcircuit_temp)
        oracle.scale(0.3)
        gpe_oracle = Group(oracle, Text("Oracle").scale(0.4)).arrange(buff=0.2, direction=DOWN)
        gpe_oracle.shift(UP*1.5+RIGHT*2)

        diffuser = Tex(r"""
        \scalebox{1.0}{
        \Qcircuit @C=1.0em @R=0.2em @!R { \\
        	 	\nghost{ {input}_{0} :  } & \lstick{ {input}_{0} :  } & \gate{\mathrm{H}} & \gate{\mathrm{X}} & \qw & \ctrl{1} & \gate{\mathrm{X}} & \gate{\mathrm{H}} & \qw & \qw & \qw\\
        	 	\nghost{ {input}_{1} :  } & \lstick{ {input}_{1} :  } & \gate{\mathrm{H}} & \gate{\mathrm{X}} & \qw & \ctrl{1} & \gate{\mathrm{X}} & \gate{\mathrm{H}} & \qw & \qw & \qw\\
        	 	\nghost{ {input}_{2} :  } & \lstick{ {input}_{2} :  } & \gate{\mathrm{H}} & \gate{\mathrm{X}} & \qw & \ctrl{1} & \gate{\mathrm{X}} & \gate{\mathrm{H}} & \qw & \qw & \qw\\
        	 	\nghost{ {input}_{3} :  } & \lstick{ {input}_{3} :  } & \gate{\mathrm{H}} & \gate{\mathrm{X}} & \qw & \ctrl{1} & \gate{\mathrm{X}} & \gate{\mathrm{H}} & \qw & \qw & \qw\\
        	 	\nghost{ {input}_{4} :  } & \lstick{ {input}_{4} :  } & \gate{\mathrm{H}} & \gate{\mathrm{X}} & \gate{\mathrm{H}} & \targ & \gate{\mathrm{H}} & \gate{\mathrm{X}} & \gate{\mathrm{H}} & \qw & \qw\\
        \\ } }
        """, tex_template=qcircuit_temp)
        diffuser.scale(0.5)
        gpe_diffuser = Group(diffuser, Text("Diffuseur").scale(0.4)).arrange(buff=0.2, direction=DOWN)
        gpe_diffuser.shift(DOWN*1.5+RIGHT*2)

        self.play(FadeIn(empty_sudoku))

        self.next_slide()

        self.play(FadeIn(gpe_oracle))

        self.next_slide()

        self.play(FadeIn(gpe_diffuser))

        self.next_slide()

        num_rep = Tex(r"Nombre de répétitions : $t = \left \lfloor \frac{\pi}{4} \sqrt{\frac{N}{m}} \right \rfloor = \left \lfloor \frac{\pi}{4} \sqrt{2^5} \right \rfloor = 4$")
        num_rep.scale(0.5)
        num_rep.shift(DOWN*0.5+LEFT*4)

        self.play(FadeIn(num_rep))

        self.next_slide()

        simulation = Text("Simulé 2048 fois : 2045 fois 01110", t2c={"2045": GREEN, "2048":RED, "01110": YELLOW})
        simulation.scale(0.3)
        simulation.shift(DOWN*1.5+LEFT*4)

        self.play(FadeIn(simulation))

        self.next_slide()

        self.play(Transform(empty_sudoku, fill_sudoku))

        self.next_slide()

        self.play(*[FadeOut(mob)for mob in self.mobjects])


class RevueOrdinateurQuantiqueSlide(Slide):
    def construct(self):
        title = Text("Ordinateur quantique / classique")
        title.shift(UP*3)

        self.play(FadeIn(title))

        self.next_slide()

        list_quant = Tex(r"\begin{itemize} \item[$\bullet$] Qubit \\ \begin{itemize} \item[$\rightarrow$] Déterministe \\ \item[$\rightarrow$] Détruit à la mesure \end{itemize} \item[$\bullet$] Portes \\ \begin{itemize} \item[$\rightarrow$] Opérateurs \\ \item[$\rightarrow$] Superposition et intrication \end{itemize} \item[$\bullet$] Algorithmes \\ \begin{itemize} \item[$\rightarrow$] Opérations parallèles \\ \item[$\rightarrow$] Résultat dans la mesure \end{itemize} \item[$\bullet$] Hardware\\ \begin{itemize} \item[$\rightarrow$] Supraconducteurs \end{itemize} \item[$\bullet$] Avantages \\ \begin{itemize} \item[$\rightarrow$] Cas spécifiques \end{itemize} \end{itemize}")
        list_quant.scale(0.4)
        list_quant.shift(LEFT*3)

        self.play(Write(list_quant))

        vertical_line = Line(UP*2.5, DOWN*3)

        self.play(Write(vertical_line))

        list_class = Tex(r"\begin{itemize} \item[$\bullet$] Bit \\ \begin{itemize} \item[$\rightarrow$] Temporaire \\ \item[$\rightarrow$] Valeur connue \end{itemize} \item[$\bullet$] Portes \\ \begin{itemize} \item[$\rightarrow$] Logique \\ \item[$\rightarrow$] Sortie unique \end{itemize} \item[$\bullet$] Algorithmes \\ \begin{itemize} \item[$\rightarrow$] Plus intuitif \\ \item[$\rightarrow$] Très généraux \end{itemize} \item[$\bullet$] Hardware\\ \begin{itemize} \item[$\rightarrow$] Transistors \end{itemize} \item[$\bullet$] Avantages \\ \begin{itemize} \item[$\rightarrow$] Cas généraux \end{itemize} \end{itemize}")
        list_class.scale(0.4)
        list_class.shift(RIGHT*3)

        self.play(Write(list_class))

        self.next_slide()

        self.play(FadeOut(title, list_quant, list_class, vertical_line))


class AvenirSlide(Slide):
    def construct(self):
        title = Text("Perspectives pour l'avenir")
        title.shift(UP*3)

        self.play(FadeIn(title))

        self.next_slide()

        list_futur = Tex(r"\begin{itemize} \item[$\bullet$] Actuel \\ \begin{itemize} \item[$\rightarrow$] Cryptographie \\ \item[$\rightarrow$] Simulation \\ \item[$\rightarrow$] ... \\ \end{itemize} \item[$\bullet$] Domaines \\ \begin{itemize} \item[$\rightarrow$] Intelligence artificiel \\ \item[$\rightarrow$] Finance \\ \item[$\rightarrow$] Pharmaceutique \\ \item[$\rightarrow$] ... \\ \end{itemize} \item[$\bullet$] Recherche \\ \begin{itemize} \item[$\rightarrow$] Qubit \\ \item[$\rightarrow$] Algorithmes \\ \item[$\rightarrow$] Réduction d'erreur \\ \item[$\rightarrow$] ... \\ \end{itemize} \item[$\bullet$] Autre \\ \begin{itemize} \item[$\rightarrow$] Éducation \\ \item[$\rightarrow$] ... \\ \end{itemize} \end{itemize}")

        list_futur.scale(0.4)
        list_futur.shift(LEFT*4+DOWN*0.5)

        self.play(Write(list_futur))

        self.next_slide()

        diagram_axes = Axes(x_range=[0, 10], y_range=[0, 10], x_length=4, y_length=5, axis_config={"include_ticks": False})
        diagram_axes.shift(2*RIGHT + DOWN*0.5)

        classical_growth = diagram_axes.plot(lambda x: 1.6**x - 1, x_range=[0, np.log(11)/np.log(1.6)], color=YELLOW)
        error_mitigation = diagram_axes.plot(lambda x: 1.3**x, x_range=[0, np.log(10)/np.log(1.3)], color=RED)
        perfect_qubit = diagram_axes.plot(lambda x: 0.2*x + 5, x_range=[0, 10], color=GREEN)

        labels_axes = diagram_axes.get_axis_labels(x_label=Tex(r"Taille\\ du calcul").scale(0.5), y_label=Tex(r"Temps\\ d'exécution").scale(0.5))

        labels_graph = Group((Tex(r"Classique").scale(0.5)).set_color(YELLOW), (Tex(r"Quantique\\ (avec erreur)").scale(0.5)).set_color(RED), (Tex(r"Quantique\\ (parfait)").scale(0.5)).set_color(GREEN)).arrange(buff=0.1, direction=DOWN)
        labels_graph.next_to(diagram_axes, RIGHT)

        actuel_dot = Dot(diagram_axes.c2p(1, 1.6**1 - 1)).set_color(PURPLE)
        arrow_to_actuel = Arrow(diagram_axes.c2p(3, 1.6**1 - 1), diagram_axes.c2p(1, 1.6**1 - 1), buff=0.1).set_color(PURPLE)
        actuel_label = (Tex(r"Actuellement").scale(0.5).next_to(arrow_to_actuel, RIGHT*0.1)).set_color(PURPLE)
        actuel_gpe = Group(actuel_dot, arrow_to_actuel, actuel_label)

        near_futur_dot = Dot(diagram_axes.c2p(4, 1.3**4)).set_color(BLUE)
        arrow_to_near_futur = Arrow(diagram_axes.c2p(6, 1.3**4), diagram_axes.c2p(4, 1.3**4), buff=0.1).set_color(BLUE)
        near_futur_label = (Tex(r"Futur proche ?").scale(0.5).next_to(arrow_to_near_futur, RIGHT*0.1)).set_color(BLUE)
        near_futur_gpe = Group(near_futur_dot, arrow_to_near_futur, near_futur_label)

        graph_gpe = Group(classical_growth, error_mitigation, perfect_qubit, labels_graph, diagram_axes, labels_axes)

        self.play(FadeIn(graph_gpe))

        self.next_slide()

        self.play(FadeIn(actuel_gpe))

        self.next_slide()

        self.play(FadeIn(near_futur_gpe))

        self.next_slide()

        self.play(FadeOut(title, list_futur, graph_gpe, actuel_gpe, near_futur_gpe))


class EndOfPresentation(Slide):
    def construct(self):
        
        remerciements = Text("Merci de votre attention !", t2c={"Merci": RED, "attention": BLUE})
        remerciements.scale(1.5)

        self.play(FadeIn(remerciements))

        self.next_slide()

        self.play(FadeOut(remerciements))

        self.next_slide()

        question = Text("Des questions ?", t2c={"questions": YELLOW})
        question.scale(1.5)

        self.play(FadeIn(question))

        self.next_slide()

        self.play(FadeOut(question))
