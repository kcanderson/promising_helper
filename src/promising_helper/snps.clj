(ns promising-helper.snps
  (:use [promising-helper.core])
  (:require [promising-helper.input :as input])
  (:require [clj-http.client :as client])  
  (:require [clojure.data.json :as json])  
  (:require [com.climate.claypoole :as claypoole]))

;; Number of allowable concurrent connections to ENSEMBL
(def NUM_ENSEMBL_CONNECTIONS 7)
(def ensembl_base "http://rest.ensembl.org/")

;; (defn- expand-region
;;   [x y left_fn right_fn]
;;   [(min (left_fn l1) (left_fn l2)) (max r1 r2)])

;; (defn- intersecting-region
;;   [regions [left right]]
;;   (first (doall (filter (fn [[l r]]
;;                           (and (<= left r)
;;                              (>= right l))) regions))))

;; (defn- expand-regions
;;   [regions]
;;   (reduce (fn [v r]
;;             (let [e (intersecting-region v r)]
;;               (if (nil? e)
;;                 (conj v r)
;;                 (conj (remove #(= e %) v)
;;                       (expand-region e r)))))
;;           []
;;           regions))

(defn gene-annotation
  [ensembl_gene_id]
  (let [url (str ensembl_base "/lookup/id/"
                 ensembl_gene_id
                 "?content-type=application/json")
        ret (json/read-str (:body (client/get url)))]
    ret))

(defn gene-location
  [ensembl_gene_id]
  (let [d (gene-annotation ensembl_gene_id)]
    (select-keys d ["seq_region_name" "start" "end" "strand"])))

(defn linked-snps
  [snp_id r_squared]
  (let [population "1000GENOMES:phase_3:KHV"
        opts (format "content-type=application/json;window_size=100;r2=%s"
                     population r_squared)
        url (format "%s/ld/human/%s/%s?%s" ensembl_base snp_id population opts)]
    ;;(try)
    (cons snp_id
          (map #(if (= snp_id (% "variation1"))
                  (% "variation2")
                  (% "variation1"))
               (json/read-str (:body (client/get url)))))
    ;;(catch Exception e nil)
    ))


(defn snp-location
  [snp_id]
  (let [url (str ensembl_base "/variation/human/"
                 snp_id
                 "?content-type=application/json")
        ret (json/read-str (:body (client/get url)))]
    (select-keys (first (ret "mappings")) ["seq_region_name" "start"])))

(defn snp-annotation
  [snp_id]
  (let [url (str ensembl_base "/variation/human/"
                 snp_id
                 "?content-type=application/json")
        ret (json/read-str (:body (client/get url)))]
    (first (ret "mappings"))))


(defn snp-annotations
  [snp_ids]
  (json/read-str
   (:body
    (client/post (str ensembl_base "variation/homo_sapiens/")
                 {:headers {"Content-Type" "application/json"
                            "Accept" "application/json"}
                  :form-params {"ids" snp_ids}
                  :content-type :json}))))

(defn ld-region
  [snp_id r_squared]
  (let [snps (remove nil? (linked-snps snp_id r_squared))
        annotations (snp-annotations snps)
        positions (remove nil?
                          (map (fn [[snp v]]
                                 (get-in v ["mappings" 0 "start"]))
                               annotations))]
    (if (not (empty? snps)) 
      {:snp (list snp_id)
       :chromosome (get-in (second (first annotations)) ["mappings" 0 "seq_region_name"])
       :start (- (reduce min positions) 50000)
       :end (+ (reduce max positions) 50000)
       })))

(defn- genes-in-region
  [chr start end]
  (let [url (str ensembl_base "/overlap/region/human/"
                 chr ":" start "-" end
                 "?feature=gene;content-type=application/json")]
    (json/read-str (:body (client/get url)))))

(defn genes-in-snp-ld-region
  [snp_id r_squared]
  (let [region (ld-region snp_id r_squared)]
    (genes-in-region (region :chromosome)
                     (region :start)
                     (region :end))))

;; (defn regions-to-genesets
;;   [regions]
;;   (into {}
;;         (mapcat (fn [[chr rgs]]
;;                   (map (fn [{start :start end :end snp :snp}]
;;                          [snp (into #{} (map #(% "external_name")
;;                                              (genes-in-region chr start end)))])
;;                        rgs))
;;                 regions)))


(defn regions-to-genesets
  [regions]
  (into {}
        (map (fn [{start :start end :end snp :snp chrom :chromosome}]
               [snp (into #{} (map #(% "external_name")
                                   (genes-in-region chrom start end)))])
             regions)))


(defn write-gmt-file
  "GMT: Gene Matrix Transposed file format.
  Tab delimited.
  First column is the name.
  Second column is a short description.
  Remaining columns are for the gene names."
  [genesets out_filename]
  (with-open [wrtr (clojure.java.io/writer out_filename)]
    (doseq [[snps gset] genesets]
      (.write wrtr (format "%s\t%s\t%s\n"
                           (first snps)
                           (str "SNPS in locus:" (clojure.string/join ", " snps))
                           (clojure.string/join "\t" gset))))))

;; (defn write-genesets-file
;;   [genesets out_filename]  
;;   (with-open [wrtr (clojure.java.io/writer out_filename)]
;;     (doseq [[snps gset] genesets]
;;       (.write wrtr (format "%s\t%s\n"
;;                            (clojure.string/join ", " snps)
;;                            (clojure.string/join ", " gset))))))

(defn- intersecting-region
  "Assumes region and other_regions are on the same chromosome."
  [region other_regions]
  (let [s1 (region :start)
        e1 (region :end)]
    (first
     (filter (fn [{s2 :start e2 :end}]
               (and (>= e1 s2) (<= s1 e2)))
             other_regions))))

(defn- expand-region
  [flank region]
  (-> region
     (assoc :start (- (region :start) flank))
     (assoc :end (+ (region :end) flank))))

(defn- expand-regions
  [flank regions]
  (map (partial expand-region flank) regions))

(defn combine-intersecting-regions
  [regions]
  (->>
   (reduce (fn [rgs r]
             (let [chrom (r :chromosome)]
               (let [sub_rgs (get rgs chrom #{})
                     int_rg (intersecting-region r sub_rgs)]
                 (assoc rgs chrom (if (nil? int_rg)
                                    (conj sub_rgs r)
                                    (conj (disj sub_rgs int_rg)
                                          (-> int_rg
                                             (assoc :snp (concat (int_rg :snp) (r :snp)))
                                             (assoc :start (min (int_rg :start) (r :start)))
                                             (assoc :end (max (int_rg :end) (r :end)))))
                                    )))))
           {}
           regions)
   (vals)
   (apply concat)))

(defn snps-to-regions 
  [snp_ids r_squared flank]
  (let [snp_ids (filter #(= "rs" (subs % 0 2)) (into #{} snp_ids))
        ;; ld_fn #(let [ret (ld-region % r_squared)]
        ;;          (println (format "%s " %))
        ;;          (Thread/sleep 1000)
        ;;          ret)
        ;; rgs (doall (map #(try (ld-region % r_squared)
        ;;                       (catch Exception e nil)
        ;;                 snp_ids))
        ;;_ (println rgs)
        regions (expand-regions flank
                                (remove nil?
                                        (claypoole/pmap NUM_ENSEMBL_CONNECTIONS
                                                        #(try (ld-region % r_squared)
                                                              (catch Exception e nil))
                                                        snp_ids)))]
    (combine-intersecting-regions regions)))

(defn snps-to-genesets
  [snp_ids r_squared flank]
  (remove (fn [[snps genes]] (empty? genes))
          (merge-genesets (regions-to-genesets (snps-to-regions snp_ids r_squared flank)))))

;; (def foo
;;   (snps-to-genesets
;;    (clojure.string/split (slurp "/home/kc/Code/Bioinformatics/reproduce_promising/snps/ad.txt") #"\n")
;;    0.5 10000))

(defn nhgri-trait-frequencies
  [gwas]
  (frequencies (map #(get % "DISEASE/TRAIT")
                    gwas)))

(defn choose-nhgri-traits
  [gwas]
  (let [vals (into {} (map-indexed vector (sort-by second > (nhgri-trait-frequencies gwas))))]
    (println (format "%s\n\nWhat traits should be used (separate by spaces): "
                     (clojure.string/join "\n"
                                          (map (fn [i]
                                                 (let [[t n] (get vals i)]
                                                   (format "[%d] %s -- Number of SNPs: %d" i t n)))
                                               (range (count vals))))))
    (let [chosen (clojure.string/split (read-line) #"\s+")
          selected (map #(first (get vals (Integer/parseInt %))) chosen)]
      (println (format "Selected: \n%s" (clojure.string/join "\n" selected)))
      selected)))

(defn distance-to-snp
  [snp ensembl_gene_id]
  (let [gloc (gene-location ensembl_gene_id)
        sloc (snp-location snp)]
    [gloc sloc]))
