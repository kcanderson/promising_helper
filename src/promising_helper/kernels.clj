(ns promising-helper.kernels
  (:require [clojure.core.matrix :as cmat])
  ;;(:require [uncomplicate.neanderthal.core :as nthal])
  ;;(:require [uncomplicate.neanderthal.native])
  )

(cmat/set-current-implementation :clatrix)
;;(cmat/set-current-implementation :nd4j)

(defn ones-matrix
  "Returns an m*n matrix of all ones."
  [m n]
  (cmat/matrix (repeat m (repeat n 1))))

(defn degree-matrix
  [mat]
  (let [[m n] (cmat/shape mat)]
    (cmat/diagonal-matrix 
     (cmat/eseq (cmat/mmul mat (ones-matrix m 1))))))

(defn laplacian-matrix
  [mat]
  (cmat/add (cmat/mul -1 mat)
            (degree-matrix mat)))

(defn normalized-laplacian-matrix
  [mat]
  (let [d (cmat/emap! #(Math/sqrt %) (degree-matrix mat))
        v (cmat/mmul (cmat/inverse d) mat d)
        [m n] (cmat/shape mat)]
    (cmat/add (cmat/identity-matrix m)
              (cmat/mul -1 v))))

(defn adjacency-mat-kernel
  [mat]
  (fn ([] mat)
    ([i1 i2] (cmat/select mat i1 i2))))

(defn kernel-shell
  [kern]
  (fn
    ([] kern)
    ([i] (cmat/get-row kern i))
    ([i1 i2]
     (if (= i1 i2) 0 (cmat/mget kern i1 i2)))))

(defn simple-inner-product-kernel
  [mat]
  (let [[m n] (cmat/shape mat)
        kern (cmat/mmul mat (cmat/transpose mat))]
    (kernel-shell kern)))

(defn commute-time-kernel
  [mat]
  (let [[m n] (cmat/shape mat)
        v (cmat/mul (/ 1 n) (ones-matrix m n))
        kern (cmat/add v (cmat/inverse
                          (cmat/add (laplacian-matrix mat)
                                    (cmat/mul -1 v))))]
    (kernel-shell kern)))

(defn von-neumann-diffusion-kernel
  [alpha mat]  
  (let [[m n] (cmat/shape mat)
        kern (cmat/inverse
              (cmat/add (cmat/identity-matrix m)
                        (cmat/mul (- alpha) mat)))]
    (kernel-shell kern)))

(defn regularized-laplacian-kernel
  [alpha mat]
  (let [[m n] (cmat/shape mat)
        kern (cmat/inverse
              (cmat/add (cmat/identity-matrix m)
                        (cmat/mul alpha (normalized-laplacian-matrix mat))))]
    (kernel-shell kern)))

