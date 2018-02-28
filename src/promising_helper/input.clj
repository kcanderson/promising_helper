(ns promising-helper.input
  (:use [promising-helper.core])
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

(defn mat-file-index-mapping
  [mat_filename]
  (with-open [rdr (clojure.java.io/reader mat_filename)]
    (let [header (clojure.string/split (first (line-seq rdr)) #"\t")]
      (into {} (map vector header (range))))
    ))

(defn mat-file-grab-entries
  "edges: collection of [first second]"
  [mat_filename edges]
  (let [index (mat-file-index-mapping mat_filename)
        rindex (into {}  (map (fn [[a b]] [b a]) index))
        item_edges (reduce (fn [m [i1 i2]]
                             (assoc m i1 (conj (get m i1 []) (index i2))))
                           {}
                           edges)
        wanted_indices (into #{} (map index (keys item_edges)))]
    (with-open [rdr (clojure.java.io/reader mat_filename)]
      (into {}
            (reduce (comp doall concat)
                    (keep-indexed (fn [i line] (if (contains? wanted_indices i)
                                                (map (fn [[i2 s]] [[(rindex i) (rindex i2)] (Float/parseFloat s)])
                                                     (select-keys (into [] (clojure.string/split (clojure.string/trim line) #"\t")) (item_edges (rindex i))))
                                                ))
                                  (rest (line-seq rdr))))))))

(defn parse-omim-phenotype-series
  [series_lines]
  (let [info_line (clojure.string/split (first series_lines) #" - ")
        f #(clojure.string/split % #"\t")
        ks (f (second series_lines))
        genes (map #(zipmap ks (f %)) (rest (rest series_lines)))
        ;;good_genes (remove #(contains? #{\? \{ \[} (first (get % "Phenotype"))) genes)
        ]
    {:disease_name (first info_line)
     :phenotypic_series_id (second info_line)
     :gene_names (map #(first (clojure.string/split (get % "Gene/Locus") #", ")) genes)
     ;;:foo (first genes)
     }))

(defn read-omim-phenotypic-series
  [lines]
  (map parse-omim-phenotype-series
       (butlast (rest (remove (comp empty? first) (partition-by (comp empty? clojure.string/trim) lines))))))

(defn gene-chromosomal-distance
  [gene_coordinates g1 g2]
  (let [[chr1 l1] ((first gene_coordinates) g1)
        [chr2 l2] ((first gene_coordinates) g2)]
    (if (= chr1 chr2)
      (java.lang.Math/abs (- l1 l2)))))

;; (defn remove-close-genes
;;   [gene_coordinates dist_threshold genes]
;;   (let [real_genes (keep #(if (nil? ((first gene_coordinates) %)) nil %) genes)]
;;     (loop [[g0 & grem] real_genes
;;            good_genes []]
;;       (if (nil? g0) good_genes
;;           (let [d (remove nil?
;;                           (map (partial gene-chromosomal-distance
;;                                         gene_coordinates g0)
;;                                grem))]
;;             (cond
;;               (empty? d) (recur grem (conj good_genes g0))
;;               (< (apply min d) dist_threshold) (recur grem good_genes)
;;               :else (recur grem (conj good_genes g0))))))))

(defn index-locations
  [sorted_gene_locations]
  (into {}
        (map-indexed
         (fn [i g]
           [(get-in g [:extra :gene_id]) (list (g :chromosome) i)])
         sorted_gene_locations)))


(defn sort-by-chromosome
  [gene_locations]
  (let [chr_map (group-by :chromosome gene_locations)]
    (reduce (fn [[cmap locs] [chr genes]]
              (let [s (vec (sort-by :start genes))]
                (list (merge cmap (index-locations s))
                      (assoc locs chr s))))
            (list {} {}) chr_map)))

;; (defn gtf-group-by-gene
;;   [gtf_parsed_lines]
;;   (vals
;;    (reduce (fn [m x]
;;              (let [k (x :ensembl_gene_id)
;;                    pid (x :protein_id)
;;                    v (get m k (-> x
;;                                  (assoc :exons [])
;;                                  (assoc :protein_ids #{})
;;                                  (dissoc :protein_id)))]
;;                (assoc m k (-> v
;;                              (assoc :start (min (v :start) (x :start)))
;;                              (assoc :end (max (v :end) (x :end)))
;;                              (assoc :protein_ids (conj (v :protein_ids) pid))
;;                              ;; TODO: Add exon information
;;                              ;;(assoc :exons (conj ))
;;                              ))))
;;            {}
;;            gtf_parsed_lines)))

(defn halve
  [n]
  (int (/ n 2)))

(defn grab-nearest-n
  [chromosome_map [chr index] num]
  (let [start (max 0 (- index (halve (inc num))))
        cmap ((second chromosome_map) chr)
        n (count cmap)
        gene_name_fn #(get-in % [:extra :gene_name])]
    (loop [i start
           s #{(gene_name_fn (get cmap index))}
           m [(get cmap index)]]
      (cond
        (or (>= i n) (= (count s) num)) m
        (contains? s (gene_name_fn (get cmap i))) (recur (inc i) s m)
        :else (let [val (get cmap i)]
                (recur (inc i) (conj s (gene_name_fn val)) (conj m val)))))))

(defn get-locus
  [chromosome_map num]
  (let [[mapper coords] chromosome_map]
    (fn [ensembl_gene_id]
      (let [loc (mapper ensembl_gene_id)]
        (if (not (nil? loc))
          (grab-nearest-n chromosome_map loc num))))))

(defn add-loci
  [gene_coordinates locus_size max_num_loci omim]
  (let [genes (take max_num_loci (shuffle (omim :ensembl_gene_ids)))
        loci (map (comp (fn [locus] (map #(get-in % [:extra :gene_id]) locus))
                        (get-locus gene_coordinates locus_size))
                  genes)]
    (assoc omim :loci loci)))

(defn omim-add-ensembl-ids-xform
  [name_ensembl_map]
  (map #(assoc % :ensembl_gene_ids
               (into #{}  (remove nil? (map name_ensembl_map (get % :gene_names)))))))

(defn omim-filter-few-genes-xform
  [min_num_genes]
  (filter #(>= (count (get % :loci)) min_num_genes)))

(defn add-loci-xform
  [gene_coords locus_size max_num_loci]
  (map (partial add-loci gene_coords locus_size max_num_loci)))

(def merge_loci_xform
  (let [f (fn [loci] (vals (merge-genesets
                           (into {}  (map-indexed (fn [i l] [[i] (into #{} l)]) loci)))))]
    (map #(assoc % :loci (f (% :loci))))
    ))

(defn omim-entries-add-loci
  [ensembl_coords ensembl_name2g_map locus_size omim_entries]
  (sequence (comp (omim-add-ensembl-ids-xform ensembl_name2g_map)
                  (add-loci-xform ensembl_coords locus_size 30)
                  merge_loci_xform
                  (omim-filter-few-genes-xform 3))
            omim_entries))

;;(println (take 1 (omim-entries-add-loci bar ensembl_name2g_map 10 omim)))

;; (with-open [rdr (clojure.java.io/reader "../../annotations/Homo_sapiens.GRCh37.70.with.entrezid.gtf")]
;;   (def foo (gtf-read-lines (line-seq rdr))))

;; (with-open [rdr (clojure.java.io/reader "../../omim/phenotypic-series-all.txt")]            
;;   (def omim (read-omim-phenotypic-series (line-seq rdr))))
;; (def ensembl_name2g_map
;;   (into {} (map (fn [x] [(get-in x [:extra :gene_name]) (get-in x [:extra :gene_id])]) foo))
;;   )
;; (def bar (sort-by-chromosome foo))
;; (sequence (comp (omim-add-ensembl-ids-xform ensembl_name2g_map)
;;                 (omim-filter-few-genes-xform 3)
;;                 (add-loci-xform bar 10 20)
;;                 merge-loci-xform)
;;           (take 10 omim))

