(ns promising-helper.main
  (:gen-class)
  (:use [promising-helper.core])
  (:require [clojure.java.io :as io])
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
  (let [protein_mapping (with-open [rdr (io/reader annotations_filename)]
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
  [interaction_filename kernel_fn]
  (with-open [rdr (clojure.java.io/reader interaction_filename)]
    (let [lines (remove (fn [[i1 i2 score]] (= i1 i2))
                        (map #(let [[i1 i2 score] (clojure.string/split % #"\t")]
                                [i1 i2 (Float/parseFloat score)])
                             (rest (line-seq rdr))))
          {mat :matrix mapping :mapping} (time (interactions-to-largest-connected-matrix lines))]
      {:kernel (kernel_fn mat)
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
      (let [kern_fn (case (options :type)
                      "reglaplacian" #(kernel/regularized-laplacian-kernel (options :alpha) %)
                      "vonneumann" #(kernel/von-neumann-diffusion-kernel (options :alpha) %)
                      "commutetime" #(kernel/commute-time-kernel %)
                      "adjacency" #(kernel/adjacency-mat-kernel %)
                      "shortestpath" #(kernel/shortest-path-kernel %)
                      "expdiffusion" #(kernel/exponential-diffusion-kernel (options :alpha) %)
                      "lapexpdiffusion" #(kernel/laplacian-exponential-diffusion-kernel (options :alpha) %)
                      )
            {kern :kernel mapping :mapping} (time (make-kernel-from-interaction-file
                                                   (options :input) kern_fn))
            name_fn (zipmap (vals mapping) (keys mapping))]
        (with-open [wrtr (clojure.java.io/writer (options :output))]
          (input/write-matrix-to-file (kern) wrtr name_fn))))))

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
   ["-c" "--common COMMON_GENES_FILE" "File containing genes to be considered in validation"]
   ["-h" "--help"]])

(defn ranked-genes-in-set
  [s ranked_list]
  (filter #(contains? s %) ranked_list))

(defn validate
  [results_file monarch_file iterations used_genes]
  (let [truth (into #{} (clojure.string/split (slurp monarch_file) #"\s+"))
        ranked_genes (evaluation/ranked-genes results_file)
        common_ranked_genes (ranked-genes-in-set (into #{} used_genes) ranked_genes)]
    {"total genes" (count ranked_genes)
     "truth genes" (count truth)
     "intersection truth" (count (clojure.set/intersection truth (into #{} ranked_genes)))
     "common genes for GSEA testing" (count common_ranked_genes)
     "common genes intersection truth" (count (clojure.set/intersection truth (into #{} common_ranked_genes)))
     "p-val" (evaluation/gsea-p-value truth common_ranked_genes iterations)
     "enrichment score" (evaluation/gsea-enrichment-score truth common_ranked_genes)}))

(defn validate-cmd
  [& args]
  (let [{:keys [options arguments errors summary]} (parse-opts args validate_opts)]
    (if (options :help)
      (println summary)
      (let [c (into #{} (clojure.string/split (slurp (options :common)) #"\s+"))
            v (validate (options :results) (options :truth) (options :pval) c)]
        (with-open [wrtr (clojure.java.io/writer (options :output))]
          (doseq [[k vv] v]
            (.write wrtr (str k "\t" vv "\n")))
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
   ["-p" "--pval PVALUE" :default 5e-8 :parse-fn #(Float/parseFloat %)]
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
                v (filter #(<= (Float/parseFloat (% "P-VALUE")) (options :pval)) v)
                snps (mapcat #(clojure.string/split (% "SNPS") #"; ") v)]
            (spit (options :output)
                  (clojure.string/join "\n" snps))
            ))))))

(defn match-results-monarch
  [results_dir monarch_dir validation_dir]
  (let [vf (fn [re path]
             (filter (fn [f] (re-find re (.getPath f)))
                     (file-seq (clojure.java.io/file path))))
        fmap (fn [files] (group-by #(first (clojure.string/split (.getName %) #"[\._]")) files))
        results (fmap (vf #".tsv" results_dir))
        monarch (fmap (vf #".txt" monarch_dir))
        validation (fmap (vf #"_common.glist" validation_dir))
        file_fn #(clojure.string/join "_" (rest (clojure.string/split % #"_")))]
    (into {} (map (fn [[k v]]
                    [k {:truth (.getPath (first v))
                        :validation (into {} (map #(vector (second (clojure.string/split (.getName %) #"_"))
                                                           (.getPath %))
                                                  (validation k)))
                        :results (into {} (map #(vector (file_fn (clojure.string/replace (.getPath %) ".tsv" ""))
                                                        (.getPath %))
                                               (results k)))}])
                  monarch))))

(def comparison_opts
  [["-r" "--results RESULTS_DIRECTORY" "directory to search for results" :default ""]
   ["-t" "--truth TRUTH_DIRECTORY" "ground truth directory"]
   ["-v" "--validation VALIDATION_DIRECTORY" "directory with common gene lists"]
   ["-o" "--output OUTPUT_FILE" "output file" :default ""]
   ["-h" "--help"]])

(defn comparison-cmd
  [& args]
  (let [{:keys [options arguments errors summary]} (parse-opts args comparison_opts)]
    (if (options :help)
      (println summary)
      (let [matches (match-results-monarch (options :results) (options :truth) (options :validation))
            ks (sort (keys (:results (second (first matches)))))]
        (with-open [wrtr (clojure.java.io/writer (options :output))]
          (.write wrtr (str (clojure.string/join "\t" (cons "disease/trait" ks)) "\n"))
          (doseq [[k {t :truth r :results v :validation}] matches]
            (let [truth (into #{} (clojure.string/split (slurp t) #"\s+"))]
              (.write wrtr (str (clojure.string/join
                                 "\t" (cons k (map #(let [r (if (get r %) (evaluation/ranked-genes (get r %)))
                                                          common_filename (get v (second (clojure.string/split % #"_")))
                                                          common (into #{} (clojure.string/split (slurp common_filename) #"\s+"))
                                                          r_common (ranked-genes-in-set common r)]
                                                      (evaluation/gsea-enrichment-score truth r_common))
                                                   ks)))
                                "\n")))))))))

(defn node-degrees-from-network
  [network_filename]
  (with-open [rdr (clojure.java.io/reader network_filename)]
    (reduce (fn [acc [g1 g2 amt]]
              (let [amt (Float/parseFloat amt)
                    v1 (get acc g1 00)
                    v2 (get acc g2 0.0)]
                (-> acc
                   (assoc g1 (+ v1 amt))
                   (assoc g2 (+ v2 amt)))))
            {}
            (map #(clojure.string/split % #"\t")
                 (rest (line-seq rdr))))))

(def degree_groups_opts
  [["-n" "--network NETWORK_FILENAME" "file name for degree groups" :default ""]
   ["-o" "--output OUTPUT_FILE" "output file" :default ""]
   ["-s" "--size DEGREE_GROUP_SIZE" "number of genes per group"
    :default 1000 :parse-fn #(Integer/parseInt %)]
   ["-h" "--help"]])

(defn degree-groups-cmd
  [& args]
  (let [{:keys [options arguments errors summary]} (parse-opts args degree_groups_opts)]
    (if (options :help)
      (println summary)
      (let [degs (node-degrees-from-network (options :network))
            sdegs (sort-by second degs)
            groups (partition (options :size) (options :size) []  sdegs)]
        (with-open [wrtr (clojure.java.io/writer (options :output))]
          (doseq [[i grp] (map-indexed vector groups)]
            (.write wrtr (clojure.string/join "\t" (concat [(str "Group " i)
                                                            (format "Contains degrees %.1f to %.1f"
                                                                    (second (first grp))
                                                                    (second (last grp)))]
                                                           (map first grp))))
            (.write wrtr "\n")))))))

(defn common-genes-among-results
  [& results_files]
  (apply clojure.set/intersection (map (comp (partial into #{}) evaluation/ranked-genes) results_files)))

(def common_opts
  [["-o" "--output OUTPUT_FILE" "output file" :default ""]
   ["-h" "--help"]])

(defn common-cmd
  [& args]
  (let [{:keys [options arguments errors summary]} (parse-opts args common_opts)]
    (if (options :help)
      (println summary)
      (try
        (spit (options :output)
              (clojure.string/join "\n" (apply common-genes-among-results arguments)))
        (catch Exception e
          (println "Error running command!\n" (.getMessage e) "\n" summary))))))

(def entries_from_kernel_opts
  [["-m" "--matrix MATRIX_FILENAME" "Kernel mat file"]
   ["-e" "--edges EDGES_FILENAME"]
   ["-o" "--output OUTPUT_FILENAME" "Output file for edges grabbed from kernel"]
   ["-h" "--help"]])

(defn entries-from-kernel-cmd
  [& args]
  (let [{:keys [options arguments errors summary]} (parse-opts args entries_from_kernel_opts)]
    (if (options :help)
      (println summary)
      (let [edges (map #(take 2 (clojure.string/split % #"\t"))
                       (clojure.string/split (slurp (options :edges)) #"\n"))
            _ (println (take 10 edges))
            entries (input/mat-file-grab-entries (options :matrix) edges)
            _ (println (take 5 entries))]
        (with-open [wrtr (clojure.java.io/writer (options :output))]
          (doseq [[[i1 i2] s] entries]
            (.write wrtr (str (clojure.string/join "\t" [i1 i2 s]) "\n"))
            ))))))

(def omim_opts
  [["-i" "--input PHENOTYPIC_SERIES" "OMIM phenotypic series file"]
   ["-n" "--num GENES_PER_LOCUS" "Number of genes per locus" :parse-fn #(Integer/parseInt %) :default 25]
   ["-a" "--annotations GTF_ANNOTATIONS" "Gene annotations"]
   ["-g" "--geneset GENESET_DIRECTORY" "Output directory to place gmt files"]
   ["-t" "--truth TRUTH_DIRECTORY" "Output directory to place true gene txt files"]
   ["-h" "--help"]])

(defn omim-cmd
  [& args]
  (let [{:keys [options arguments errors summary]} (parse-opts args omim_opts)]
    (if (options :help)
      (println summary)
      (let [ensembl (with-open [rdr (io/reader (options :annotations))]
                      (input/gtf-read-lines (line-seq rdr)))
            ensembl_coords (input/sort-by-chromosome ensembl)
            omim (with-open [rdr (clojure.java.io/reader (options :input))]            
                   (input/read-omim-phenotypic-series (line-seq rdr)))
            ensembl_name2g_map (into {} (map (fn [x] [(get-in x [:extra :gene_name])
                                                     (get-in x [:extra :gene_id])])
                                             ensembl))
            loci_omim (input/omim-entries-add-loci ensembl_coords ensembl_name2g_map (options :num) omim)
            g2n (into {} (map (fn [[a b]] [b a]) ensembl_name2g_map))]
        (doseq [o loci_omim]
          (let [filename (-> (str "omim-" (o :disease_name))
                            (clojure.string/replace "(" "-")
                            (clojure.string/replace ")" "-")
                            (clojure.string/replace "/" "-")
                            (clojure.string/replace " " "-")
                            )]
            (spit (str (options :truth) (java.io.File/separator) filename ".txt")
                  (clojure.string/join "\n" (o :gene_names)))
            (with-open [wrtr (io/writer (str (options :geneset) (java.io.File/separator) filename ".gmt"))]
              (doseq [[i l] (map-indexed vector (o :loci))]
                (.write wrtr (clojure.string/join "\t" (concat [(str (o :disease_name) " locus " i)
                                                           (str "Locus for "
                                                                (clojure.string/join ", " (map g2n (clojure.set/intersection l (into #{} (o :ensembl_gene_ids))))))]
                                                               (remove empty? (map g2n l)))))
                (.write wrtr "\n")))))))))
        
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
                                                           ["enrichment-figure" "Make enrichment plots from results"]
                                                           ["comparison" "Compare all methods"]
                                                           ["commonalities" "Find common genes among results"]
                                                           ["degree-groups" "Make degree groups GMT"]
                                                           ["entries-from-kernel" "Grab edges from kernel"]
                                                           ["omim-genesets-cmd" "Make genesets from OMIM phenotypic series"]])]
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
      :comparison (apply comparison-cmd args)
      :commonalities (apply common-cmd args)
      :degree-groups (apply degree-groups-cmd args)
      :entries-from-kernel (apply entries-from-kernel-cmd args)
      :omim-genesets-cmd (apply omim-cmd args)
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





