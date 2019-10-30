# from: http://www.gnuplotting.org/code/default_color_map3.gnu
unset key
set cbtics scale 0
set palette defined ( 0 "#000090",\
                      1 "#000fff",\
                      2 "#0090ff",\
                      3 "#0fffee",\
                      4 "#90ff70",\
                      5 "#ffee00",\
                      6 "#ff7000",\
                      7 "#ee0000",\
                      8 "#7f0000")

plot "/tmp/image.txt" with image, "/tmp/evaluations.txt" pt 7 ps 1, "/tmp/evaluations.txt" with lines ls 0.1
pause -1
