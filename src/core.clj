(ns core
  (:require [clojure.core.matrix :as m]
            [utils :as utils]
            [types :refer [->Node ->Element ->Condition]]))

(m/set-current-implementation :vectorz)

(defn run [opts]
  (let [mesh (utils/read-input "./data/viga_15_prueba_y_error.dat")
        params (-> mesh :params)
        nodes (-> mesh :nodes)
        elements (-> mesh :elements)
        conditions (-> mesh :conditions)
        local-k-array (map #(utils/create-local-k % (-> params :EI)) elements)
        local-b-array (map #(utils/create-local-b % (-> params :f)) elements)
        {K_temp :K b_temp :b} (utils/assembly nodes elements local-k-array local-b-array)
        b_neumann (utils/apply-neumann-condition (-> conditions :neumann) b_temp)
        {K :K b :b} (utils/apply-dirichlet-condition (-> conditions :dirichlet :x) K_temp b_neumann)
        X (utils/calculate-fem K b)
        final-X (utils/insert-condition-values X (-> conditions :dirichlet :x))]
    (m/pm final-X)
    (utils/write-to-file "./data/viga_15_prueba_y_error.post.res" final-X)))
