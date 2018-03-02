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

(defn matrix-pow
  [mat p]
  (iterate ))

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
                        (cmat/mul alpha (laplacian-matrix mat))))]
    (kernel-shell kern)))

(defn exponential-diffusion-kernel
  [alpha mat]
  (let [kern (matrix-exp (cmat/mul alpha mat))]
    (fn
      ([] kern)
      ([i1 i2] (cmat/select kern i1 i2)))))

(defn laplacian-exponential-diffusion-kernel
  [alpha mat]
  (let [kern (matrix-exp (cmat/mul (- alpha)
                                   (laplacian-matrix mat)))]
    (fn [i1 i2]
      (cmat/select kern i1 i2))))

(defn cosine-kernel
  [mat]
  (let [[m n] (cmat/shape mat)]
    (cmat/matrix
     (partition
      n (for [i (range m)
              j (range n)]
          (let [a (cmat/get-row mat i)
                b (cmat/get-row mat j)]
            (cmat/div
             (cmat/dot a b)
             (java.lang.Math/sqrt
              (* (cmat/dot a a) (cmat/dot b b))))))))))

(defn queue
  [coll]
  (into clojure.lang.PersistentQueue/EMPTY coll))

(def infinity Integer/MAX_VALUE)

(defn update-neighbors
  [q dists curr neighbors]
  (let [new_dist (inc (get dists curr))]
    (reduce (fn [[q' dists'] i]
              (if (< new_dist (dists i))
                [(conj q' i) (assoc dists' i new_dist)]
                [q' dists']))
            [q dists]
            neighbors)))

(defn get-neighbors
  [mat index]
  (keep (fn [x] (if (> (second x) 0) (first x)))
        (map-indexed list
                     (cmat/eseq (cmat/get-row mat index)))))

(defn propogate-distances-from-start
  [mat]
  (let [[m n] (cmat/shape mat)
        neighbor_fn (memoize (partial get-neighbors mat))]
    (fn [kern start]
      (cmat/set-selection! kern start start 0)
      (loop [q (queue (list start))]
        (cond
          (empty? q) kern
          :else (let [curr (peek q)
                      neighbors (neighbor_fn curr)]
                  (recur
                   (reduce
                    (fn [q' n]
                      (let [curr_dist (cmat/select kern start n)
                            new_dist (+ (cmat/select kern start curr)
                                        (cmat/select mat curr n))]
                        (if (< new_dist curr_dist)
                          (do (cmat/set-selection! kern start n new_dist)
                              ;;(cmat/set-selection! kern n start new_dist)
                              (conj q' n))
                          q')))
                    (pop q) neighbors))))))))

(defn shortest-path-kernel
  ([mat] (shortest-path-kernel #(/ 1.0 (+ 1.0 %)) mat))
  ([dist_to_sim_fn mat]
   (let [[m n] (cmat/shape mat)
         start_fn (propogate-distances-from-start mat)
         kern (reduce (fn [k start]
                        (do (println start)
                            (start_fn k start)))
                      (cmat/add infinity (cmat/zero-matrix m n))
                      (range m))
         kern (cmat/emap! dist_to_sim_fn kern)]
     (fn
       ([] kern)
       ([i1 i2] (cmat/select kern i1 i2))))))

(defn p-step-random-walk-kernel
  [p a mat]
  (let [[m n] (cmat/shape mat)
        b (cmat/add (cmat/mul a (cmat/identity-matrix m))
                    (cmat/mul -1.0 (normalized-laplacian-matrix mat)))]
    
    ))
