(ns promising-paper.input
  (:use [promising-paper.core])
  (:require [clojure.core.matrix :as cmat]))

(defn parse-tsv
  [lines]
  (let [f #(clojure.string/split % #"\t")
        ks (f (first lines))]
    (map #(zipmap ks (f %)) (rest lines))))


;; Read string network
(defn interactome-line
  [line]
  (let [[p1 p2 score] (clojure.string/split line #" ")
        f #(subs % 5)]
    (list (list (f p1) (f p2)) (/ (read-string score) 1000.0))))

(defn combined-score
  [line]
  (let [[p1 p2 & vals] (clojure.string/split line #" ")
        f #(subs % 5)]
    (list (list (f p1) (f p2))
          (/ (float (Integer/parseInt (last vals))) 1000.0))))

(defn noisy-or
  [coll]
  (- 1 (reduce * (map #(- 1 %) coll))))

(defn denoisy-or
  [total denoisy_amount]
  (+ 1.0 (/ (- total 1.0) (- 1.0 denoisy_amount))))

(defn interactome-parse-line
  [line]
  (let [[p1 p2 neighborhood fusion
         cooccurrence coexpression
         experimental database textmining combined] (clojure.string/split line #" ")
        f #(subs % 5)]
    (list (list (f p1) (f p2))
          (zipmap [:neighborhood :fusion :cooccurrence :coexpression
                   :experimental :database :textmining :combined]
                  (map (comp #(/ % 1000.0) #(Integer/parseInt %))
                       (list neighborhood fusion cooccurrence coexpression
                             experimental database textmining combined))))
    ))

(defn sans-lit-line
  [line]
  (let [[p1 p2 neighborhood fusion
         cooccurrence coexpression
         experimental database textmining combined] (clojure.string/split line #" ")
        f #(subs % 5)]
    (list (list (f p1) (f p2))
          (noisy-or
           (map (comp #(/ % 1000.0) float #(Integer/parseInt %))
                (list neighborhood fusion cooccurrence coexpression experimental))))))

(defn sans-comention-line
  [line]
  (let [[p1 p2 neighborhood fusion
         cooccurence coexpression
         experimental database textmining combined] (clojure.string/split line #" ")
        f #(subs % 5)]
    (list (list (f p1) (f p2))
          (noisy-or
           (map (comp #(/ % 1000.0) read-string)
                (list neighborhood fusion cooccurence coexpression experimental database))))))

(defn experimental-only-line
  [line]
  (let [[p1 p2 neighborhood fusion
         cooccurence coexpression
         experimental database textmining combined] (clojure.string/split line #" ")
        f #(subs % 5)]
    (list (list (f p1) (f p2))
          (/ (read-string experimental) 1000.0))))

(defn coexpression-only-line
  [line]
  (let [[p1 p2 neighborhood fusion
         cooccurence coexpression
         experimental database textmining combined] (clojure.string/split line #" ")
        f #(subs % 5)]
    (list (list (f p1) (f p2))
          (/ (read-string coexpression) 1000.0))))

(defn database-and-textmining-line
  [line]
  (let [[p1 p2 neighborhood fusion
         cooccurence coexpression
         experimental database textmining combined] (clojure.string/split line #" ")
        f #(subs % 5)]
    (list (list (f p1) (f p2))
          (noisy-or [(/ (read-string textmining) 1000.0)
                     (/ (read-string database) 1000.0)]))))

(defn textmining-only-line
  [line]
  (let [[p1 p2 neighborhood fusion
         cooccurence coexpression
         experimental database textmining combined] (clojure.string/split line #" ")
        f #(subs % 5)]
    (list (list (f p1) (f p2))
          (/ (read-string textmining) 1000.0))))

(defn database-only-line
  [line]
  (let [[p1 p2 neighborhood fusion
         cooccurence coexpression
         experimental database textmining combined] (clojure.string/split line #" ")
        f #(subs % 5)]
    (list (list (f p1) (f p2))
          (/ (read-string textmining) 1000.0))))

(defn cooccurrence-only-line
  [line]
  (let [[p1 p2 neighborhood fusion
         cooccurence coexpression
         experimental database textmining combined] (clojure.string/split line #" ")
        f #(subs % 5)]
    (list (list (f p1) (f p2))
          (/ (read-string cooccurence) 1000.0))))

(defn fusion-only-line
  [line]
  (let [[p1 p2 neighborhood fusion
         cooccurence coexpression
         experimental database textmining combined] (clojure.string/split line #" ")
        f #(subs % 5)]
    (list (list (f p1) (f p2))
          (/ (read-string fusion) 1000.0))))

(defn no-neighborhood-or-lit-line
  [line]
  (let [[p1 p2 neighborhood fusion
         cooccurence coexpression
         experimental database textmining combined] (clojure.string/split line #" ")
        f #(subs % 5)]
    (list (list (f p1) (f p2))
          (noisy-or
           (map (comp #(/ % 1000.0) read-string)
                (list fusion cooccurence coexpression experimental))))))

(defn no-specific-comentions-line
  [ensembl_protein_ids]
  (let [s (into #{} ensembl_protein_ids)]
    (fn [line]
      (let [[p1 p2 neighborhood fusion
             cooccurence coexpression
             experimental database textmining combined] (clojure.string/split line #" ")
            f #(subs % 5)]
        (if (and (contains? s (f p1)) (contains? s (f p2)))
          (sans-comention-line line)
          (list (list (f p1) (f p2)) (/ (Double/parseDouble combined) 1000.0)))))))

(defn neighborhood-only-line
  [line]
  (let [[p1 p2 neighborhood fusion
         cooccurence coexpression
         experimental database textmining combined] (clojure.string/split line #" ")
        f #(subs % 5)]
    (list (list (f p1) (f p2))
          (noisy-or
           (map (comp #(/ % 1000.0) read-string)
                (list neighborhood))))))

(defn read-interactome
  ([lines threshold]
   (read-interactome lines threshold interactome-line))
  ([lines threshold score_fn]
   (let [interactions (map score_fn lines)]
     (filter (comp #(> % threshold) second)
             interactions))))

(defn populate-network
  [interactions protein_mapper]
  (let [n (count protein_mapper)
        mat (cmat/zero-matrix n n)]
    (doseq [[[p1 p2] score] interactions]
      (let [v1 (protein_mapper p1)
            v2 (protein_mapper p2)]
        (cmat/mset! mat v1 v2 score)
        (cmat/mset! mat v2 v1 score)))
    mat))

(defn protein-set
  [interactions]
  (reduce
   (fn [s [[i j] score]]
     (-> s (conj i) (conj j)))
   #{} interactions))

(defn interactions-to-network
  [protein_interactions]
  (let [pset (time (protein-set protein_interactions))
        pmapper (zipmap pset (range))]
    {:protein_map pmapper
     :matrix (time (populate-network protein_interactions pmapper))}))

(defn read-interaction-file
  ([filename score_threshold] (read-interaction-file filename score_threshold interactome-line ))
  ([filename score_threshold score_fn]
   (with-open [rdr (clojure.java.io/reader filename)]
     (interactions-to-network (read-interactome (rest (line-seq rdr)) score_threshold score_fn))
     )))

(defn read-matrix
  [lines]
  (let [f #(clojure.string/split % #"\t")
        genes (f (first lines))]
    {:matrix (cmat/matrix
              (map (fn [line]
                     (map #(Float/parseFloat %) (f line)))
                   (rest lines)))
     :genename_map (into {} (map-indexed (fn [i g] [g i]) genes))}))



(defn legitimate-chromosome-name
  [name]
  (< (count name) 3))

(defn gtf-get-chromosome
  [chr_str]
  (let [s (first (clojure.string/split chr_str #"_"))]
    (if (= (apply str (take 5 s)) "HSCHR")
      (subs s 5)
      s)))

(defn gtf-parse-extra-fields
  [extra]
  (let [s (clojure.string/split extra #";\s*")
        pairs (map #(let [[a b] (clojure.string/split % #" ")]
                      [(keyword a) (subs b 1 (dec (count b)))]) s)]
    (into {} pairs)))

(defn gtf-line
  [split_line]
  (let [extra (clojure.string/trim (nth split_line 8))]
    {:chromosome (gtf-get-chromosome (first split_line)) ;;(first split_line)
     :start (Integer/parseInt (nth split_line 3))
     :end (Integer/parseInt (nth split_line 4))
     :extra (gtf-parse-extra-fields extra)
     }))

(defn gtf-read-lines
  [gtf_lines]
  (transduce (comp (map #(clojure.string/split % #"\t"))
                   (filter #(and (= (nth % 2) "CDS")
                               (legitimate-chromosome-name (first %))))
                   (map gtf-line))
             conj
             (drop-while #(= \# (first %)) gtf_lines)))

(defn write-matrix-to-file
  ([mat writer] (write-matrix-to-file mat writer nil))
  ([mat writer name_fn]
   (let [[m n] (cmat/shape mat)]
     (if-not (nil? name_fn)
       (.write writer (str (clojure.string/join "\t" (map name_fn (range m)))
                           "\n")))
     (doseq [i (range n)]
       (let [row (assoc (into [] (cmat/get-row mat i)) i 0)]
         (.write writer (str (clojure.string/join "\t" (map (comp str float) row))
                             "\n")))))))

