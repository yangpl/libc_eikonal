../bin/stencil_test  <<EOD
1
100 0              ! point A
0 100              ! point B
0 0                ! point C
200 200              ! source
1000 500 0 0         ! vitesse A
1000 500 0 0         ! vitesse B
1000 500 0 0         ! vitesse C
-1                 ! norder
EOD

../bin/stencil_test  <<EOD
1
100 0              ! point A
0 100              ! point B
0 0                ! point C
-200 200              ! source
1000 500 0 0         ! vitesse A
1000 500 0 0         ! vitesse B
1000 500 0 0         ! vitesse C
-1                 ! norder
EOD

../bin/stencil_test  <<EOD
1
100 0              ! point A
0 100              ! point B
0 0                ! point C
-200 -200              ! source
1000 500 0 0         ! vitesse A
1000 500 0 0         ! vitesse B
1000 500 0 0         ! vitesse C
-1                 ! norder
EOD

../bin/stencil_test  <<EOD
1
100 0              ! point A
0 100              ! point B
0 0                ! point C
200 -200              ! source
1000 500 0 0         ! vitesse A
1000 500 0 0         ! vitesse B
1000 500 0 0         ! vitesse C
-1                 ! norder
EOD

