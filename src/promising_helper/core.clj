(ns promising-helper.core
  (:require [clojure.core.matrix :as cmat]))


(defn populate-matrix
  [interactions mapping]
  (let [n (count mapping)
        mat (cmat/zero-matrix n n)]
    (doseq [[i1 i2 score] interactions]
      (let [v1 (mapping i1)
            v2 (mapping i2)]
        (cmat/mset! mat v1 v2 score)
        (cmat/mset! mat v2 v1 score)))
    mat))

(defn- unique-items
  [interactions]
  (reduce
   (fn [s [i j score]]
     (-> s (conj i) (conj j)))
   #{} interactions))

(defn interactions-to-matrix
  "Takes an iterable of the form: [item1 item2 score]"
  [interactions]
  (let [items (unique-items interactions)
        mapping (zipmap items (range))]
    {:mapping mapping
     :matrix (populate-matrix interactions mapping)}))



(defn get-neighbors
  [mat index]
  (keep (fn [x] (if (> (second x) 0) (first x)))
        (map-indexed list
                     (cmat/eseq (cmat/get-row mat index)))))

(defn connected-component
  ([mat starting_index]
   (connected-component mat starting_index (partial get-neighbors mat)))
  ([mat starting_index neighbor_fn]
   (loop [nodes #{}
          [next & rem] (list starting_index)]
     (cond
       (nil? next) nodes
       :else (let [neighbors (remove #(contains? nodes %)
                                     (neighbor_fn next))]
               (recur (conj nodes next) (into () (into #{} (concat rem neighbors)))))))))

(defn connected-components
  [mat]
  (let [[m n] (cmat/shape mat)
        neighbor_fn (memoize (partial get-neighbors mat))]
    (loop [left (into #{} (range m))
           components (list)]
      (cond
        (empty? left) components
        :else (let [seed (first left)
                    component (connected-component mat seed neighbor_fn)]
                (recur (clojure.set/difference left component)
                       (conj components component)))))))

(defn largest-connected-component-of-matrix
  [mat]
  (let [main_component (apply max-key count (connected-components mat))
        indices (into [] (sort main_component))
        ret (cmat/zero-matrix (count indices) (count indices))]
    {:matrix (do (doseq [i (range (count indices))
                         j (range (count indices))]
                   (cmat/mset! ret i j (cmat/mget mat (indices i) (indices j))))
                 ret)
     :mapping (zipmap indices (range))}))

(defn connected-components-from-interactions
  [interactions]
  (let [f (fn [comps x] (first (filter #(contains? % x) comps)))]
    (loop [[[curr_i curr_j curr_score] & rem] interactions
           components #{}]
      (if (nil? curr_i) components
          (let [connected_i (f components curr_i)
                connected_j (f components curr_j)]
            (cond
              (and (nil? connected_i) (nil? connected_j)) (recur rem (conj components #{curr_i curr_j}))
              (nil? connected_i) (recur rem (-> (disj components connected_j)
                                               (conj (conj connected_j curr_i))))
              (nil? connected_j) (recur rem (-> (disj components connected_i)
                                               (conj (conj connected_i curr_j))))
              (= connected_i connected_j) (recur rem components)
              :else (recur rem  (-> components
                                   (disj connected_i)
                                   (disj connected_j)
                                   (conj (clojure.set/union connected_i connected_j))))
              ))))))

(defn largest-connected-component-from-interactions
  [interactions]
  (apply max-key count (connected-components-from-interactions interactions)))

(defn interactions-to-largest-connected-matrix
  [interactions]
  (let [component (largest-connected-component-from-interactions interactions)
        mapping (zipmap component (range))]
    {:matrix (populate-matrix (filter (fn [[i1 i2 score]]
                                        (or (contains? component i1)
                                           (contains? component i2)))
                                      interactions)
                              mapping)
     :mapping mapping}))

