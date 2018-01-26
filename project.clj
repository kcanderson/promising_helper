(defproject promising_helper "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :main promising-helper.main
  :aot :all
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.8.0"]
                 [com.climate/claypoole "1.1.4"]                 
                 [org.clojure/data.json "0.2.6"]
                 [incanter "1.9.0"]
                 [clj-http "3.7.0"]
                 [cheshire "5.8.0"]
                 [org.clojure/tools.cli "0.3.5"]
                 [clj-sub-command "0.3.0"]
                 [net.mikera/core.matrix "0.61.0"]
                 ;;[cav/mtj "0.4.2-KC-PATCH"]
                 [clatrix "0.5.0"]
                 ;;[nd4clj "0.1.0-SNAPSHOT"]
                 ;;[uncomplicate/neanderthal "0.18.0"]
                 
                 ]
  :jvm-opts ["-Xmx16g" "-server"])








