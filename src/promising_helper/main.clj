(ns promising-helper.main
  (:gen-class)
  (:use [promising-helper.core])
  (:require [promising-helper.snps :as snp])
  (:require [promising-helper.input :as input])
  (:require [promising-helper.kernels :as kernel])
  (:require [promising-helper.evaluation :as evaluation])
  (:require [clj-sub-command.core :refer [sub-command candidate-message]])
  (:require [clojure.tools.cli :refer [parse-opts]]))

(defn third [coll]
  (second (rest coll)))

(defn make-genesets
  ([snps_filename out_filename]
   (make-genesets snps_filename
                  out_filename
                  0.5))
  ([snps_filename out_filename r_squared]
   (make-genesets snps_filename
                  out_filename
                  r_squared
                  10000))
  ([snps_filename out_filename r_squared flank]
   (println (format "---- Reading SNPs from %s ----" snps_filename))
   (let [snps (clojure.string/split (slurp snps_filename) #"\s+")]
     (println "---- Generating genesets for each SNP. This will query ENSEMBL servers for 1000 Genomes LD data ----")
     (let [gsets (snp/snps-to-genesets snps r_squared flank)]
       (println (format "---- Writing genesets to %s ----" out_filename))
       (snp/write-gmt-file gsets out_filename)))))

(def gsets_opts
  [["-f" "--flank FLANK" "Flank size in bases" :default 10000 :parse-fn #(Integer/parseInt %)]
   ["-r" "--rsquared R2" "r^2 for SNPs in 1000 genomes" :default 0.5 :parse-fn #(Float/parseFloat %)]
   ["-h" "--help"]])

(defn genesets-cmd
  [& args]
  (let [{:keys [options arguments errors summary]} (parse-opts args gsets_opts)]
    (if (options :help)
      (println summary)
      (make-genesets (first arguments)
                     (second arguments)
                     (options :rsquared)
                     (options :flank)))))

(defn make-string-derived-network
  [string_filename threshold include_textmining? annotations_filename out_filename]
  (println (format "------ Reading annotations from %s ------" annotations_filename))
  (let [protein_mapping (with-open [rdr (clojure.java.io/reader annotations_filename)]
                          (into {} (map (fn [x] [(get-in x [:extra :protein_id])
                                                (get-in x [:extra :gene_name])])
                                        (input/gtf-read-lines (line-seq rdr)))))
        line_fn (if include_textmining? input/combined-score input/sans-lit-line)]
    (with-open [rdr (clojure.java.io/reader string_filename)
                wrtr (clojure.java.io/writer out_filename)]
      (.write wrtr (format "gene1\tgene2\tscore\n"))
      (doseq [[[p1 p2] score] (input/read-interactome (rest (line-seq rdr)) threshold line_fn)]
        (.write wrtr (format "%s\t%s\t%f\n" (protein_mapping p1) (protein_mapping p2) score))))))

(defn make-pf-derived-network
  [pf_filename annotations_filename out_filename]
  (println (format "------ Reading annotations from %s ------" annotations_filename))
  (let [mapping (with-open [rdr (clojure.java.io/reader annotations_filename)]
                  (into {} (map (fn [x] [(get-in x [:extra :entrez_id])
                                        (get-in x [:extra :gene_name])])
                                (remove #(nil? (get-in % [:extra :entrez_id]))
                                        (input/gtf-read-lines (line-seq rdr))))))
        f (comp mapping #(subs % 11 (dec (count %))))]
    (with-open [rdr (clojure.java.io/reader pf_filename)
                wrtr (clojure.java.io/writer out_filename)]
      (.write wrtr (format "gene1\tgene2\tscore\n"))
      (doseq [line (rest (line-seq rdr))]
        (let [[a b] (clojure.string/split line #",")]
          (if (and (not (nil? (f a))) (not (nil? (f b))))
            (.write wrtr (format "%s\t%s\t1.0\n" (f a) (f b)))))))))

(def network_opts
  [["-t" "--type NETWORK_TYPE" "network type: string, pf, or " :default "string" :parse-fn identity]
   ["-a" "--annotations GTF_FILE" "GTF ENSEMBL annotations" :default "" :parse-fn identity]
   ["-m" "--textmining" "Include text mining" :default false]
   ["-i" "--input INPUT_NETWORK" "Input network file" :default "" :parse-fn identity]
   ["-o" "--output OUTPUT_NETWORK" "File name for output network" :default "" :parse-fn identity]
   ["-s" "--threshold THRESHOLD" "Minimum score for network" :default 0.15 :parse-fn #(Float/parseFloat %)]
   ["-h" "--help"]])

(defn derive-network-cmd
  [& args]
  (let [{:keys [options arguments errors summary]} (parse-opts args network_opts)]
    (if (options :help)
      (println summary)
      (case (options :type)
        "string" (make-string-derived-network (options :input)
                                              (options :threshold)
                                              (options :textmining)
                                              (options :annotations)
                                              (options :output))
        "pf" (make-pf-derived-network (options :input)
                                      (options :annotations)
                                      (options :output))))))


(defn make-kernel-from-interaction-file
  [interaction_filename alpha]
  (with-open [rdr (clojure.java.io/reader interaction_filename)]
    (let [lines (map #(let [[i1 i2 score] (clojure.string/split % #"\t")]
                        [i1 i2 (Float/parseFloat score)])
                     (rest (line-seq rdr)))
          {mat :matrix mapping :mapping} (time (interactions-to-largest-connected-matrix lines))]
      {:kernel (kernel/regularized-laplacian-kernel alpha mat)
       :mapping mapping})))

(def kernel_opts
  [["-t" "--type KERNEL" "network type: reglaplacian " :default "reglaplacian"]
   ["-i" "--input INPUT_NETWORK" "Input network file in tsv format" :default ""]
   ["-o" "--output OUPUT_MATRIX" "File name for kernel in matrix format" :default ""]
   ["-a" "--alpha ALPHA_PARAMETER" "Alpha parameter for some kernels" :default 0.001 :parse-fn #(Float/parseFloat %)]
   ["-h" "--help"]])

(defn kernel-cmd
  [& args]
  (let [{:keys [options arguments errors summary]} (parse-opts args kernel_opts)]
    (if (options :help)
      (println summary)
      (case (options :type)
        "reglaplacian" (let [{kern :kernel mapping :mapping} (time (make-kernel-from-interaction-file
                                                                    (options :input)
                                                                    (options :alpha)))
                             name_fn (zipmap (vals mapping) (keys mapping))]
                         (with-open [wrtr (clojure.java.io/writer (options :output))]
                           (input/write-matrix-to-file (kern) wrtr name_fn)))))))

;; MONARCH
(def monarch_opts
  [["-i" "--input MONARCH_FILE" "Input monarch tab-separated value file" :default ""]
   ["-o" "--output OUTPUT_FILE" "Output filename" :default ""]
   ["-h" "--help"]])

(defn monarch-cmd
  [& args]
  (let [{:keys [options arguments errors summary]} (parse-opts args monarch_opts)]
    (if (options :help)
      (println summary)
      (spit (options :output)
            (clojure.string/join "\n" (evaluation/monarch-causal-genes (options :input)))))))

;; (def validate_opts
;;   [["-r" "--results RESULTS_DIR" "Directory containg prioritization results" :default ""]
;;    ["-t" "--truth GROUND_TRUTH_DIR" "Directory containing known true genes files" :default ""]
;;    ["-o" "--output " "File name for output" :default ""]
;;    ["-p" "--pval ITERATIONS" "Permutations to run for p-val" :parse-fn #(Integer/parseInt %) :default 10000]
;;    ["-h" "--help"]])

(def validate_opts
  [["-r" "--results RESULTS_FILE" "File with results"]
   ["-t" "--truth GROUND_TRUTH" "File with ground truth genes"]
   ["-o" "--output OUTPUT_FILE" "File for output"]
   ["-p" "--pval ITERATIONS" "Number of permutations for p-value" :default 1000 :parse-fn #(Integer/parseInt %)]
   ["-h" "--help"]])

(defn match-results-monarch
  [results_dir monarch_dir]
  (let [vf (fn [re path]
             (filter (fn [f] (re-find re (.getPath f)))
                     (file-seq (clojure.java.io/file path))))
        fmap (fn [files] (group-by #(first (clojure.string/split (.getName %) #"\.")) files))
        results (fmap (vf #".tsv" results_dir))
        monarch (fmap (vf #".txt" monarch_dir))]
    (into {} (map (fn [[k v]]
                    [k [(.getPath (first v)) (map #(.getPath %) (results k))]])
                  monarch))))

(defn validate
  [results_file monarch_file iterations]
  (let [truth (into #{} (clojure.string/split (slurp monarch_file) #"\s+"))
        ranked_genes (evaluation/ranked-genes results_file)]
    (println (count truth) (count ranked_genes) (count (clojure.set/intersection truth (into #{} ranked_genes))))
    (println (take 10 ranked_genes))
    {"total genes" (count ranked_genes)
     "truth genes" (count truth)
     "intersection truth" (count (clojure.set/intersection truth (into #{} ranked_genes)))
     "p-val" (evaluation/gsea-p-value truth ranked_genes iterations)
     "enrichment score" (evaluation/gsea-enrichment-score truth ranked_genes)}))

(defn validate-cmd
  [& args]
  (let [{:keys [options arguments errors summary]} (parse-opts args validate_opts)]
    (if (options :help)
      (println summary)
      (let [v (validate (options :results) (options :truth) (options :pval))]
        (with-open [wrtr (clojure.java.io/writer (options :output))]
          (.write wrtr (str "total genes\t" (v "total genes") "\n"))
          (.write wrtr (str "num truth genes\t" (v "truth genes") "\n"))
          (.write wrtr (str "num intersection\t" (v "intersection truth") "\n"))
          (.write wrtr (str "enrichment score\t" (v "enrichment score") "\n"))
          (.write wrtr (str "p-val\t" (v "p-val") "\n"))
          )))))

(def gseaplot_opts
  [["-r" "--results RESULTS_FILE" "File with results"]
   ["-t" "--truth GROUND_TRUTH" "File with ground truth genes"]
   ["-o" "--output OUTPUT_PNG_FILE" "File for output"]
   ["-h" "--help"]])

(defn make-enrichment-plot
  [truth_filename all_result_filenames output_filename]
  (let [truth (into #{} (clojure.string/split (slurp truth_filename) #"\s+"))
        ;;f #(apply str (rest (clojure.string/split % (re-pattern (java.io.File/separator)))))
        f #(clojure.string/replace % (re-pattern (str ".*results" (java.io.File/separator))) "")
        r (into {} (map #(vector (f %) (evaluation/ranked-genes %)) all_result_filenames))]
    (evaluation/make-enrichment-figure truth r output_filename)))

(defn enrichment-cmd
  [& args]
  (let [{:keys [options arguments errors summary]} (parse-opts args gseaplot_opts)]
    (if (options :help)
      (println summary)
      (make-enrichment-plot (options :truth) arguments (options :output)))))

;; (defn validate-cmd
;;   [& args]
;;   (let [{:keys [options arguments errors summary]} (parse-opts args validate_opts)
;;         type_fn (fn [path]
;;                   (let [v (clojure.string/split path (re-pattern (java.io.File/separator)))]
;;                     (clojure.string/join ", " (take-last 2 (butlast v)))))
;;         validate_results_fn (fn [tr all_rs]
;;                               (into {} (map (fn [r]
;;                                               [(type_fn r) (validate r tr (options :pval))])
;;                                             all_rs)))]
;;     (if (options :help)
;;       (println summary)
;;       (with-open [wrtr (clojure.java.io/writer (options :output))]
;;         (let [matches (match-results-monarch (options :results) (options :truth))
;;               ks (map type_fn (second (second (first matches))))
;;               header (cons "disease/phenotype" ks)]
;;           (.write wrtr (str (clojure.string/join "\t" header) "\n"))
;;           (doseq [[k [t rs]] matches]
;;             (let [rvalidates (validate_results_fn t rs)]
;;               (.write wrtr (str (clojure.string/join "\t" (cons k (map rvalidates ks))) "\n")))))))
;;     (System/exit 0)))

;; SNPs
(def select_traits_opts
  [["-i" "--input SNPS_FILE" :default ""]
   ["-o" "--output OUTPUT_TRAIT_FILE" :default ""]
   ["-h" "--help"]])

(defn select-traits-cmd
  [& args]
  (let [{:keys [options arguments errors summary]} (parse-opts args select_traits_opts)]
    (if (options :help)
      (println summary)
      (with-open [rdr (clojure.java.io/reader (options :input))]
        (spit (options :output)
              (clojure.string/join "\n" (snp/choose-nhgri-traits (input/parse-tsv (line-seq rdr)))))))))

(def snps_opts
  [["-i" "--input SNPS_FILE" :default ""]
   ["-t" "--traits TRAITS_FILE" :default ""]
   ["-o" "--output OUTPUT_SNP_FILE" :default ""]
   ["-h" "--help"]])

(defn snps-cmd
  [& args]
  (let [{:keys [options arguments errors summary]} (parse-opts args snps_opts)]
    (if (options :help)
      (println summary)
      (let [traits (into #{} (clojure.string/split (slurp (options :traits)) #"\n"))]
        (println (slurp (options :traits)))
        (with-open [rdr (clojure.java.io/reader (options :input))]
          (let [v (filter #(contains? traits (% "DISEASE/TRAIT"))
                          (input/parse-tsv (line-seq rdr)))
                v (filter #(not= "NR" (% "SNPS")) v)
                snps (mapcat #(clojure.string/split (% "SNPS") #"; ") v)]
            (spit (options :output)
                  (clojure.string/join "\n" snps))
            ))))))

(defn -main
  "Testing, 1, 2, 3!"
  [& args]
  (let [[opts cmd args help cands] (sub-command args
                                                "Usage yo"
                                                :options [["-h" "--help" "Show help" :default false :flag true]]
                                                :commands [["genesets" "Make gene sets from SNPs list"]
                                                           ["network" "Make network"]
                                                           ["kernel" "Kernelize interaction network"]
                                                           ["monarch" "Parse MONARCH data to grab causal relations"]
                                                           ["validate" "Perform GSEA on results"]
                                                           ["select-traits" "Select important traits from NHGRI GWAS table"]
                                                           ["snps" "Pull out desired SNPs from NHGRI GWAS table"]
                                                           ["enrichment-figure" "Make enrichment plots from results"]])]
    (when (:help opts)
      (println help)
      (System/exit 0))
    (case cmd
      :genesets (apply genesets-cmd args)
      :network (apply derive-network-cmd args)
      :kernel (apply kernel-cmd args)
      :monarch (apply monarch-cmd args)
      :validate (apply validate-cmd args)
      :select-traits (apply select-traits-cmd args)
      :snps (apply snps-cmd args)
      :enrichment-figure (apply enrichment-cmd args)
      (println (str "Invalid command. See 'foo --help'.\n\n"
                    (candidate-message cands))))
    (System/exit 0)))

;; (def foo
;;   (with-open [rdr (clojure.java.io/reader "/home/kc/Code/Bioinformatics/interactome/musings/data/ncbi_ids.txt")]
;;     (into {}
;;           (map #(let [x (clojure.string/split % #"\t")]
;;                   [(second x) (try (Integer/parseInt (last x))
;;                                    (catch Exception e nil))])
;;                (rest (line-seq rdr))))))

;; (with-open [rdr (clojure.java.io/reader "/home/kc/Code/Bioinformatics/interactome/musings/data/Ensembl/Homo_sapiens.GRCh37.70.gtf")
;;             wrtr (clojure.java.io/writer "/home/kc/Code/Bioinformatics/interactome/musings/data/Ensembl/Homo_sapiens.GRCh37.70.with.entrezid.gtf")]
;;   (doseq [l (map (fn [line]
;;                    (let [s (re-find #"gene_name \S+;" line)
;;                          s (subs s 0 (dec (count s)))
;;                          [a b] (clojure.string/split s #" ")
;;                          gene (subs b 1 (dec (count b)))]
;;                      (str line " entrez_id \"" (foo gene) "\";\n")))
;;                  (line-seq rdr))]
;;     (.write wrtr l)))


;; (with-open [rdr (clojure.java.io/reader "/home/kc/Code/Bioinformatics/reproduce_promising/annotations/../annotations/Homo_sapiens.GRCh37.70.with.hgnc.gtf")]
;;   (let [foo (input/gtf-read-lines (take 100000 (line-seq rdr)))]
;;     (println (first foo)))
;;   ;;(input/gtf-read-lines )
;;   )





