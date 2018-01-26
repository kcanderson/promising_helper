(ns promising-paper.evaluation
  (:require [promising-paper.input :as input])
  ;;(:require [incanter.charts :as charts])
  ;;(:require [incanter.core])
  )

(defn positions
  [pred coll]
  (keep-indexed (fn [idx x]
                  (when (pred x) idx))
                coll))

(defn ranked-genes
  [results_filename]
  (with-open [rdr (clojure.java.io/reader results_filename)]
    (doall
     (map #(get % "gene") (input/parse-tsv (line-seq rdr))))))

(defn gsea-curve
  [truth ranked_list]
  (let [truth_set (into #{} truth)
        intsct (clojure.set/intersection truth_set (into #{} ranked_list))
        delta_up (/ 1.0 (count intsct))
        delta_down (- (/ 1.0 (- (count ranked_list) (count intsct))))]
    (reduce (fn [acc g]
              (conj acc (+ (last acc)
                           (if (contains? truth_set g)
                             delta_up
                             delta_down))))
            [0.0]
            ranked_list)))

(defn gsea-enrichment-score
  [truth ranked_list]
  (reduce max (gsea-curve truth ranked_list)))

(defn gsea-p-value
  [truth ranked_list num_permutations]
  (let [real (gsea-enrichment-score truth ranked_list)
        n (count ranked_list)
        iteration_fn (fn [] (reduce max (gsea-curve truth (shuffle ranked_list))))]
    (/ (float (count (filter #(< real %)
                             (repeatedly num_permutations iteration_fn))))
       num_permutations)))

(def good_relation_labels #{"pathogenic_for_condition" "contributes to" "has phenotype"})

(defn monarch-causal-genes
  [monarch_filename]
  (with-open [rdr (clojure.java.io/reader monarch_filename)]
    (into #{}
          (map #(get % "subject_label")
               (filter #(contains? good_relation_labels (get % "relation_label"))
                       (input/parse-tsv (line-seq rdr)))))))

;; (with-open [rdr (clojure.java.io/reader "../../monarch/t1d.tsv")]
;;   (def t1d_genes (into #{}
;;                        (map #(get % "subject_label")
;;                             (filter #(contains? good_relation_labels (get % "relation_label"))
;;                                     (doall (parse-tsv (line-seq rdr))))))))

;; (with-open [rdr (clojure.java.io/reader "../../monarch/t1d.tsv")]
;;   (frequencies
;;       (map #(get % "relation_label")
;;            (doall (parse-tsv (line-seq rdr))))
;;       ))

;; (println (clojure.string/join "\n" ad_genes))
;; (count (clojure.set/intersection (into #{} ranked_genes) ad_genes))

;; (gsea-enrichment-score ad_genes ranked_genes)
;; (gsea-p-value ad_genes ranked_genes )

;; (defn mean [coll] (/ (float (reduce + coll)) (count coll)))

;; (mean
;;  (map #(/ (float %) (count ranked_genes))
;;       (positions #(contains? ad_genes %) ranked_genes)))

;; (def ranked_genes_pf (ranked-genes "../../results/pf/string_notm/ad.tsv"))
;; (def ranked_genes_pr (ranked-genes "../../results/promising/string_notm/ad.tsv"))

;; (incanter.core/view
;;  (charts/add-lines
;;   (charts/xy-plot (map #(/ (float %) (count ranked_genes_pr)) (range (count ranked_genes_pr)))
;;                   (gsea-curve ad_genes ranked_genes_pr))
;;   (map #(/ (float %) (count ranked_genes_pf)) (range (count ranked_genes_pf)))
;;   (gsea-curve ad_genes ranked_genes_pf)))

;; (gsea-p-value ad_genes ranked_genes_pr 1000)
;; (gsea-p-value ad_genes ranked_genes_pf 10000)


;; (with-open [rdr (clojure.java.io/reader "/home/kc/Downloads/gwas-association-downloaded_2018-01-24-Type.tsv")]
;;   (println
;;    (clojure.string/join "\no"
;;                         (sort-by second >
;;                                  (frequencies
;;                                   (map #(get % "DISEASE/TRAIT")
;;                                        (parse-tsv (line-seq rdr))))))))

;; (with-open [rdr (clojure.java.io/reader "/home/kc/Code/Bioinformatics/reproduce_promising/source_snps/gwas-association-downloaded_2018-01-24-Type.tsv")]
;;   (spit "/home/kc/Downloads/t1d.txt"
;;    (clojure.string/join "\n"
;;                         (map #(get % "SNPS")
;;                              (filter (fn [x] (contains? #{"Type 1 diabetes"} (x "DISEASE/TRAIT")))
;;                                      (parse-tsv (line-seq rdr)))))))


;; (with-open [rdr (clojure.java.io/reader "../../genesets/t1d.gmt")]
;;   (def foo (into #{} (mapcat #(rest (rest (clojure.string/split % #"\t")))
;;                              (line-seq rdr)))))


