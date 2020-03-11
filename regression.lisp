(unlock-package :sb-ext)
(ql:quickload '(:inferior-shell :cl-csv :vgplot :inferior-shell))

(defpackage :corona
  (:use :cl :inferior-shell :cl-csv :alexandria
   :vgplot))

(in-package :corona)

(load "polynomial-regression.lisp")
;; --------------------------------------------------------------------------------
;; Coronavirus interpretation
;; --------------------------------------------------------------------------------
(defparameter corona (make-hash-table :test 'equal))

(defun confirmed-case (row)
  (setf (gethash (concatenate 'string (car row) " " (cadr row)) corona)
        (let ((hs (make-hash-table)))
          (loop for x in #1=(nthcdr 4 row)
                and i from 1 to (length #1#)
                do (if (equal x "")  (gethash (- i 1) hs)
                       (setf (gethash i hs)
                             (coerce (parse-integer x) 'double-float))))
          hs)))

(defun inactive-case (row)
  (let ((hs (gethash (concatenate 'string (car row) " " (cadr row)) corona)))
    (loop for x in #1=(nthcdr 4 row)
          and i from 1 to (length #1#)
          do (setf (gethash i hs)
                   (if (equal x "")  (gethash i hs)
                       (- (gethash i hs)
                           (coerce (parse-integer x) 'double-float)))))))

;; We factor out china for two reasons:
;; Ground zero occured in centeral china. They progressed to the later stages
;; of an epidemic far before the rest of the world started to see a majority of cases.
;; China also has had one of the most effiecient lock-down we've seen, so the spread
;; has mostly been contained.
(defparameter non-chinese-cases (make-hash-table))

(defun daily-cases ()
  "daily aggregate"
  (map nil #'(lambda (x)
             (loop for z in (hash-table-alist (cdr x))
                   do (setf (gethash (car z) non-chinese-cases)
                            ;; so either there is no value currently, in which case we use 0
                            (+ (or (gethash (car z) non-chinese-cases) 0)
                                ;; or
                                (or (cdr z) 0))))) ;; or use the previous day...

       ;; (mapcar #'(lambda (x)
       ;;            (hash-table-alist (cdr x))
       ;;            )
       ;; We filter out china
       (remove-if #'(lambda (x) (cl-ppcre:scan "China" (car x)))
                  (hash-table-alist corona))))

(defun rate (data)
  "rate of change"
  (if (consp (cddr data))
      (cons (- (/ (cadr data) (car data)) 1)
            (rate (cdr data)))
      (cons (- (/ (cadr data) (car data)) 1) nil)))

(defun array-to-list (l)
  (loop for x below (car (array-dimensions l))
        collect (loop for a below (cadr (array-dimensions l))
                      collect (aref l x a))))

(defun list-to-array (lst)
  (make-array `(1 ,(length lst))
              :element-type 'double-float
              :initial-contents (list (mapcar #'(lambda (x)
                                               (coerce x 'double-float))
                                           lst))))
 
(defun polynomial-regression (data coeff)
  "Regression based on coeff."
  (let ((x-values (list-to-array (loop for x from 1 to (length data)
                                     collect (coerce x 'double-float)))))
    (alexandria:flatten (array-to-list (polynomial-fit x-values (list-to-array data) coeff)))))

(defun polynomial-create (l)
  "Creates a closure from a polynomial regression."
  (lambda (x)
    (declare (ignorable x))
    (reduce #'+ (loop for y from 0 to (- (length l) 1) and z in l
                      collect (* z (expt x y))))))

(defun extrapolate (current rate)
  "this should be our sigmoid function actually"
  ;; currently we don't know the value of L
  ;; but the curve should look like
  (* current (expt 2.718 rate)))


  ;; Setup data

(defun main ()
  (if (directory "data")
      (run/s "cd data; git pull")
      (run/s "git clone https://github.com/CSSEGISandData/COVID-19 data"))

  ;; we want current active cases. To get that we have to subtract the number of recovered + deaths from the number of conformed cases.
  ;; We start with confirmed
  (let ((cases (directory "data/csse_covid_19_data/csse_covid_19_time_series/*.csv")))
    (with-open-file (s (print (directory-namestring (find-if #'(lambda (x) (cl-ppcre:scan "Confirmed" x)) cases))))
      (map nil #'confirmed-case (cdr (read-csv s))))

  ;; subtract deaths
    (with-open-file (s (find-if #'(lambda (x) (cl-ppcre:scan "Deaths" x)) cases))
    (map nil #'inactive-case (cdr (read-csv s))))

  ;; subtract recovered
    (with-open-file (s (find-if #'(lambda (x) (cl-ppcre:scan "Recovered" x)) cases))
    (map nil #'inactive-case (cdr (read-csv s)))))

  ;; linear regression of rate of change
  (daily-cases)

  (let ((slope-regression
        (polynomial-create
         (polynomial-regression
          (rate (hash-table-values non-chinese-cases)) 1)))
      ;; we'll extrapolate out the rate of change using it's linear regression
      (extrapolated-hs (make-hash-table)))

;;  Most recent rate
    (let* ((recent-day (apply #'max (hash-table-keys non-chinese-cases)))
           (var (gethash recent-day non-chinese-cases))
           (rate (lastcar (rate(hash-table-values non-chinese-cases)))))
      ;; We'll extrapolate out 3 months
      (dotimes (x 90)
        (let ((d (extrapolate var rate)))
          (setf (gethash (+ x recent-day) extrapolated-hs) d)
          (setf var d
                rate (funcall slope-regression (+ x recent-day))))))

    (plot
        ;; current cases
     (hash-table-keys non-chinese-cases)
     (hash-table-values non-chinese-cases)

      ;; extrapolation
      (reverse (hash-table-keys extrapolated-hs))
      (reverse (hash-table-values extrapolated-hs)))))

  ;; (sb-ext:save-lisp-and-die #P"corona" :toplevel #'corona::main :executable t :compression t)

  ;; rates of change
  ;; (hash-table-keys non-chinese-cases)
  ;; (rate (hash-table-values non-chinese-cases))

;; linear regression of that rate
;; (loop for x from 1 to 60 collect x)
;; (loop for x from 1 to 60 collect (funcall rate-regression x))


;; (defun logistic-growth (current rate time)
;;   ;; capacity is the carrying capacity, the limiting value of N
;;   ;; the inflection point occurs as N = K/2
;;   ;; (x = * (+ 1 n) * p) ;; we know that the rate, p, is something like 1.15
;;   ;; (print (* current (expt 2.718 (* growth time)))))))
;;   ;; b is determined by K / N(0) - 1
;;   (let ((b (derivative (/ (- capacity current) current))))))
;;   ;; In this case we actually are trying to find capacity
;;   ;; I think....
;;   (derivative K (* current (+ 1 (* b (expt 2.718 (* (- rate) time))))))
;;   ;; I think if we can run gradient descent over which interval of t = time
;;   ;; has the steepest increase, then we can start to see when it will hit 0
;;   ;; (/ capacity (+ 1 (* b (expt 2.718 (* (- rate) time)))))
;;   ;; )
;; (gradi (rate (hash-table-values non-chinese-cases)))

