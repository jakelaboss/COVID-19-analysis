(in-package :corona)
;; --------------------------------------------------------------------------------
;; Cholesky Decomposition
;; --------------------------------------------------------------------------------

(defun cholesky! (A b &key (eps single-float-epsilon))
  "given
   [1] A (required)
       ==> vector of length N*(N+1)/2 representing
           a two dimensional N by N positive definite matrix
           with elements stored columnwise; i.e, elements with indices
           0     1     2     3     4     5    ...
           correspond to matrix positions
           1,1   1,2   2,2   1,3   2,3   3,3  ...
   [2] b (required)
       <=> right-hand vector on input;
           solution X on output
   [3] eps (keyword; single-float-epsilon)
       number used to test for a computationally singular matrix
solves A X = b using the Cholesky decomposition method, and
returns
   [1] X, the solution vector (stored in b)"
  (let* ((N (length b))
         (N-1 (1- N))
         (XL (make-array `(,N ,N)))
         (Y (make-array N))
         (D (make-array N))
         (L 0))
    ;(declare (dynamic-extent XL Y D))
    (assert (= (length A) (/ (* N (1+ N)) 2)))
    ;;; Factor into triangular and diagonal form ...
    (setf (aref D 0) (realpart (aref A 0)))
    (dotimes (i N-1)
      (dotimes (j (1+ i))
        (incf L)
        (setf (aref XL (1+ i) j) (/ (conjugate (aref A L)) (aref D j)))
        (when (not (zerop j))
          (dotimes (k j)
            (decf (aref XL (1+ i) j)
                  (/ (* (aref XL (1+ i) k)
                        (conjugate (aref XL j k))
                        (aref D k))
                     (aref D j))))))
      (incf L)
      (setf (aref D (1+ I)) (realpart (aref A L)))
      (dotimes (k (1+ I))
        (decf (aref D (1+ I))
              (* (aref D k) (expt (abs (aref XL (1+ i) k)) 2))))
      ;;; Test for nonpositive value (i.e., matrix is too close to
      ;;; being singular)
      (if (< (aref D (1+ I)) eps)
        (error "(aref D ~D) = ~F < ~F = esp~&         matrix is computationally singular"
               (1+ I) (aref D (1+ I)) eps)))
    ;;; Solve for intermediate column vector solution ...
    (setf (aref Y 0) (aref b 0))
    (dotimes (k N-1)
      (setf (aref Y (1+ k)) (aref b (1+ k)))
      (dotimes (j (1+ k))
        (decf (aref Y (1+ k))
              (* (aref XL (1+ k) j) (aref Y j)))))
    ;;; Solve for final column vector solution ...
    (setf (aref b N-1) (/ (aref Y N-1) (aref D N-1)))
    (let ((k-index N-1))
      (dotimes (k-count N-1 (values b))
        (decf k-index)
        (setf (aref b k-index) (/ (aref Y k-index) (aref D k-index)))
        (dotimes (j (- N k-index 1))
          (decf (aref b k-index)
                (* (conjugate (aref XL (+ k-index j 1) k-index))
                   (aref b (+ k-index j 1)))))))))

;; Calculates the Cholesky decomposition matrix L
;; for a positive-definite, symmetric nxn matrix A.
(defun chol (A)
  (let* ((n (car (array-dimensions A)))
         (L (make-array `(,n ,n) :element-type 'double-float :initial-element 0d0)))
    (do ((k 0 (incf k))) ((> k (- n 1)) nil)
        ;; First, calculate diagonal elements L_kk.
        (setf (aref L k k)
              (sqrt (- (aref A k k)
                       (do* ((j 0 (incf j))
                             (sum (expt (aref L k j) 2)
                                  (incf sum (expt (aref L k j) 2))))
                            ((> j (- k 1)) sum)))))
        ;; Then, all elements below a diagonal element, L_ik, i=k+1..n.
        (do ((i (+ k 1) (incf i)))
            ((> i (- n 1)) nil)
            (setf (aref L i k)
                  (/ (- (aref A i k)
                        (do* ((j 0 (incf j))
                              (sum (* (aref L i j) (aref L k j))
                                   (incf sum (* (aref L i j) (aref L k j)))))
                             ((> j (- k 1)) sum)))
                     (aref L k k)))))
    ;; Return the calculated matrix L.
    L))

;; case of matrix stored as a list of lists (inner lists are rows of matrix)
;; as above, returns the Cholesky decomposition matrix of a square positive-definite, symmetric matrix
(defun cholesky (m)
  (let ((l (list (list (sqrt (caar m))))) x (j 0) i)
    (dolist (cm (cdr m) (mapcar #'(lambda (x) (nconc x (make-list (- (length m) (length x)) :initial-element 0))) l))
      (setq x (list (/ (car cm) (caar l))) i 0)
      (dolist (cl (cdr l)) 
        (setf (cdr (last x)) (list (/ (- (elt cm (incf i)) (*v x cl)) (car (last cl))))))
      (setf (cdr (last l)) (list (nconc x (list (sqrt (- (elt cm (incf j)) (*v x x))))))))))

;; where *v is the scalar product defined as
(defun *v (v1 v2) (reduce #'+ (mapcar #'* v1 v2)))


;; --------------------------------------------------------------------------------
;; Matrix Transposition
;; --------------------------------------------------------------------------------

;; if initial data is a list we can transpose with the list function
(defun transpose (m)
  (apply #'mapcar #'list m))

;; If the matrix A is given as a 2D array:
;; Transpose a mxn matrix A to a nxm matrix B=A'.
(defun mtp (A)
  (let* ((m (array-dimension A 0))
         (n (array-dimension A 1))
         (B (make-array `(,n ,m) :element-type 'double-float :initial-element 0d0)))
    (loop for i from 0 below m do
      (loop for j from 0 below n do
        (setf (aref B j i)
              (aref A i j))))
    B))

;; --------------------------------------------------------------------------------
;; Matrix Multiplication
;; --------------------------------------------------------------------------------

(defun matrix-multiply (a b)
  (flet ((col (mat i) (mapcar #'(lambda (row) (elt row i)) mat))
         (row (mat i) (elt mat i)))
    (loop for row from 0 below (length a)
          collect (loop for col from 0 below (length (row b 0))
                        collect (apply #'+ (mapcar #'* (row a row) (col b col)))))))

;; example use:
;; (matrix-multiply '((1 2) (3 4)) '((-3 -8 3) (-2 1 4)))

(defun matrix-multiply (m1 m2)
  (mapcar
   (lambda (row)
     (apply #'mapcar
            (lambda (&rest column)
              (apply #'+ (mapcar #'* row column)))
            m2))
   m1))

;; The following version uses 2D arrays as inputs.

(defun mmul (A B)
  (let* ((m (car (array-dimensions A)))
         (n (cadr (array-dimensions A)))
         (l (cadr (array-dimensions B)))
         (C (make-array `(,m ,l) :element-type 'double-float :initial-element 0d0)))
    (loop for i from 0 to (- m 1) do
      (loop for k from 0 to (- l 1) do
        (setf (aref C i k)
              (loop for j from 0 to (- n 1)
                    sum (* (aref A i j)
                         (aref B j k))))))
    C))

;; Another version
(defun mmult (a b)
  (loop
    with m = (array-dimension a 0)
    with n = (array-dimension a 1)
    with l = (array-dimension b 1)
    with c = (make-array (list m l) :initial-element 0)
    for i below m do
      (loop for k below l do
        (setf (aref c i k)
              (loop for j below n
                    sum (* (aref a i j)
                           (aref b j k)))))
    finally (return c)))


;; --------------------------------------------------------------------------------
;; Multiple Regression
;; --------------------------------------------------------------------------------

;; Solve a linear system AX=B where A is symmetric and positive definite, so it can be Cholesky decomposed.
(defun linsys (A B)
  (let* ((n (car  (array-dimensions A)))
         (m (cadr (array-dimensions B)))
         (y (make-array n        :element-type 'double-float :initial-element 0.0d0))
         (X (make-array `(,n ,m) :element-type 'double-float :initial-element 0.0d0))
         (L (chol A))) ; A=LL'
    (loop for col from 0 to (- m 1) do
       ;; Forward substitution: y = L\B
       (loop for k from 0 to (- n 1)
             do (setf (aref y k)
                      (/ (- (aref B k col)
                            (loop for j from 0 to (- k 1)
                                  sum (* (aref L k j)
                                         (aref y j))))
                         (aref L k k))))
       ;; Back substitution. x=L'\y
       (loop for k from (- n 1) downto 0
             do (setf (aref X k col)
                      (/ (- (aref y k)
                            (loop for j from (+ k 1) to (- n 1)
                                  sum (* (aref L j k)
                                         (aref X j col))))
                         (aref L k k)))))
    X))

;; Solve a linear least squares problem. Ax=b, with A being mxn, with m>n.
;; Solves the linear system A'Ax=A'b.
(defun lsqr (A b)
  (linsys (mmul (mtp A) A)
          (mmul (mtp A) b)))

;; --------------------------------------------------------------------------------
;;; Polynomial Regression
;; --------------------------------------------------------------------------------

;; Least square fit of a polynomial of order n the x-y-curve
(defun polynomial-fit (x y n)
  (let* ((m (cadr (array-dimensions x)))
         (A (make-array `(,m ,(+ n 1)) :element-type 'double-float :initial-element 0d0)))
    (loop for i from 0 to (- m 1) do
      (loop for j from 0 to n do
        (setf (aref A i j)
              (expt (aref x 0 i) j))))
    (lsqr A (mtp y))))
