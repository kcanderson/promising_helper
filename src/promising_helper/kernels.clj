(ns promising-helper.kernels
  (:require [clojure.core.matrix :as cmat])
  ;;(:require [uncomplicate.neanderthal.core :as nthal])
  ;;(:require [uncomplicate.neanderthal.native])
  )

(cmat/set-current-implementation :mtj)
;;(cmat/set-current-implementation :clatrix)
;;(cmat/set-current-implementation :nd4j)

(defn ones-matrix
  "Returns an m*n matrix of all ones."
  [m n]
  (cmat/matrix (repeat m (repeat n 1))))

(defn my-diagonal-matrix
  [diag_vals]
  (let [m (cmat/identity-matrix (count (into [] diag_vals)))]
    (doseq [[i d] (map-indexed list diag_vals)]
      (cmat/mset! m i i d))
    m))

(defn degree-matrix
  [mat]
  (let [[m n] (cmat/shape mat)]
    (my-diagonal-matrix (cmat/get-column (cmat/mmul mat (ones-matrix m 1)) 0))))

(defn matrix-pow
  [mat p]
  (loop [n 1
         acc mat]
    (if (>= n p) acc
        (recur (inc n) (cmat/mmul acc mat)))))

(defn matrix-exp
  [mat]
  (reduce (fn [mr [m k i]]
            (do (println i)
                (cmat/add mr (cmat/div m k))))
          (apply cmat/identity-matrix (cmat/shape mat))
          (take 8
                (iterate (fn [[m k i]]
                           (let [knext (* k (inc i))]
                             (list (cmat/mmul m mat)
                                   knext
                                   (inc i))))
                         (list mat 1.0 1)))))

(defn laplacian-matrix
  [mat]
  (cmat/sub! (degree-matrix mat) mat))

(defn normalized-laplacian-matrix
  [mat]
  (let [[m n] (cmat/shape mat)
        d (let [v (degree-matrix mat)]
            (doseq [i (range m)]
              (cmat/mset! v i i (/ 1.0 (Math/sqrt (cmat/mget v i i)))))
            v)]
    (cmat/sub! (cmat/mutable (cmat/identity-matrix m))
               (cmat/mmul d mat d))))

;; (normalized-laplacian-matrix (cmat/matrix [[0 1 3]
;;                                            [1 0 1]
;;                                            [3 1 0]]))
;; (normalized-laplacian-matrix
;;  (cmat/matrix [[0 1 3]
;;                [1 0 1]
;;                [3 1 0]]))

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
                        (cmat/mul alpha (laplacian-matrix mat))))]
    (kernel-shell kern)))

(defn exponential-diffusion-kernel
  [alpha mat]
  (let [kern (matrix-exp (cmat/mul alpha mat))]
    
    (kernel-shell kern)))

(defn laplacian-exponential-diffusion-kernel
  [alpha mat]
  (let [kern (matrix-exp (cmat/mul (- alpha)
                                   (laplacian-matrix mat)))]
    (kernel-shell kern)))

;; (defn cosine-kernel
;;   [mat]
;;   (let [[m n] (cmat/shape mat)
;;         square_fn (memoize #(cmat/dot % %))]
;;     (cmat/matrix
;;      (partition
;;       n (for [i (range m)
;;               j (range n)]
;;           (let [a (cmat/get-row mat i)
;;                 b (cmat/get-row mat j)]
;;             (if (zero? j) (println i))
;;             (cmat/div
;;              (cmat/dot a b)
;;              (java.lang.Math/sqrt
;;               (* (cmat/dot a a) (cmat/dot b b)
;;                ;;(square_fn a) (square_fn b)
;;                )))))))))

(defn cosine-kernel
  [mat]
  (let [v (cmat/mmul mat (cmat/transpose mat))
        squares (cmat/emap! #(Math/sqrt %) (cmat/diagonal v))
        vv (cmat/outer-product squares squares)]
    (kernel-shell (cmat/div! v vv))))

(defn queue
  [coll]
  (into clojure.lang.PersistentQueue/EMPTY coll))

(def infinity Integer/MAX_VALUE)

(defn get-neighbors
  [mat index]
  (keep (fn [x] (if (> (second x) 0) (first x)))
        (map-indexed list
                     (cmat/eseq (cmat/get-row mat index)))))

(defn update-neighbors
  [q dists dist_fn neighbors]
  (let [curr (peek q)
        dist_to_curr (get dists curr)]
    (loop [[n0 & nrest] neighbors
           q' (pop q)]
      (if (nil? n0) q'
          (let [curr_dist (get dists n0)
                new_dist (+ dist_to_curr (dist_fn curr n0))]
            (if (< new_dist curr_dist)
              (recur nrest (do (assoc! dists n0 new_dist)
                               (conj q' n0)))
              (recur nrest q')))))))

(defn shortest-paths-from
  [num_nodes dist_fn neighbor_fn start]
  (let [dists (assoc! (transient (vec (repeat num_nodes Double/POSITIVE_INFINITY)))
                      start 0)]
    (loop [q (queue (list start))]
      (cond
        (empty? q) (persistent! dists)
        :else (recur (update-neighbors q dists dist_fn (neighbor_fn (peek q))))))))

(defn shortest-path-kernel
  ([mat] (shortest-path-kernel mat (memoize (partial get-neighbors mat))))
  ([mat neighbor_fn]
   (let [[m n] (cmat/shape mat)
         kern (cmat/mul -1 (cmat/matrix
                            (pmap (fn [i]
                                    (do (if (zero? (mod i 100)) (println i))
                                        (shortest-paths-from
                                         m
                                         (partial cmat/mget mat)
                                         neighbor_fn i)))
                                  (range m))))]
     (kernel-shell kern))))

(defn p-step-random-walk-kernel
  [p b mat]
  (let [[m n] (cmat/shape mat)
        v (cmat/sub! (cmat/mul b (cmat/identity-matrix m))
                     (normalized-laplacian-matrix mat))]
    (kernel-shell (matrix-pow v p))))

