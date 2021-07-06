(ns utils
  (:require [clojure.core.match :refer [match]]
            [clojure.core.matrix :as m]
            [clojure.math.numeric-tower :as math]
            [types :refer [->Node ->Element ->Condition]]))

(m/set-current-implementation :vectorz)

(defn get-node-by-id
  "Gets a node from <nodes> by its <id>"
  [id nodes]
  (some #(and (= id (:id %)) %) nodes))

(defn node-index "get node from index" [node]
  (- (-> node :id) 1))

(defn get-nodes-by-element "Calculate number of nodes by element"
  [element]
  (-> element keys count))

(defn get-element-by-id
  "Gets an element from <elements> by its id"
  [id elements]
  (some #(and (= id (:id %)) %) elements))

(defn add-row
  "Add row at <index> in <matrix>"
  [matrix index row]
  (let [[rows-len cols-le] (m/shape matrix)]
    (match [index]
      [0] (m/matrix (m/join-along 0 [row] matrix))
      [rows-len] (m/matrix (m/join-along 0 matrix [row]))
      :else (let [pre (m/submatrix matrix [[0 index]])
                  post (m/submatrix matrix [[index (- rows-len index)]])]
              (m/matrix (m/join-along 0 pre [row] post))))))

; (def test-mat (m/matrix [[1]
;                          [2]
;                          [3]]))
; (m/pm (add-row test-mat 3 [0]))

(defn remove-row
  "Removes row <index> from matrix"
  [matrix index]
  (let [[rows-len cols-len] (m/shape matrix)
        new-len (- rows-len 1)]
    (match [index]
      [0] (m/submatrix matrix [[1 new-len]])
      [new-len] (m/submatrix matrix [[0 new-len]])
      :else (let [pre (m/submatrix matrix [[0 index]])
                  post (m/submatrix matrix [[(inc index) (- rows-len index 1)]])]
              (m/matrix (m/join pre post))))))

(defn remove-col
  "Removes row <index> from matrix"
  [matrix index]
  (let [[rows-len cols-len] (m/shape matrix)
        new-len (- cols-len 1)]
    (match [index]
      [0] (m/submatrix matrix [nil [1 new-len]])
      [new-len] (m/submatrix matrix [nil [0 new-len]])
      :else (let [pre (m/submatrix matrix [nil [0 index]])
                  post (m/submatrix matrix [nil [(inc index) (- cols-len index 1)]])]
              (m/matrix (m/join-along 1 pre post))))))

(defn jacobian
  "Return jacobian for trans"
  [element]
  (let [{n1 :n1 n2 :n2 n3 :n3 n4 :n4} element
        {x1 :x y1 :y z1 :z} n1
        {x2 :x y2 :y z2 :z} n2
        {x3 :x y3 :y z3 :z} n3
        {x4 :x y4 :y z4 :z} n4]
    (m/det (m/matrix [[(- x2 x1) (- x3 x1) (- x4 x1)]
                      [(- y2 y1) (- y3 y1) (- y4 y1)]
                      [(- z2 z1) (- z3 z1) (- z4 z1)]]))))

(defn insert-submatrix-in-matrix
  "Insert <submatrix> at <start> in <matrix>"
  [submatrix pos matrix]
  (let [[rows-len cols-len] (m/shape submatrix)
        {abs-i :i abs-j :j} pos
        temp (m/matrix matrix)]
    (do
      (doseq [i (range rows-len)
              j (range cols-len)]
        (m/mset! temp (+ abs-i i) (+ abs-j j) (m/mget submatrix i j)))
      temp)))

(defn replace-div-by-zero "replace divide by zero" [a b]
  (if (zero? b)
    0.000001
    (/ a b)))

(defn create-local-k
  "Create Local K system"
  [element EI]
  (let [{n1 :n1 n2 :n2 n3 :n3 n4 :n4} element
        jacob (jacobian element)
        x1 (get-in element [:n1 :x])
        x2 (get-in element [:n2 :x])
        x8 (get-in element [:n8 :x])
        c1 (replace-div-by-zero 1 (math/expt (- x2 x1) 2))
        c2 (* (replace-div-by-zero 1 (- x2 x1)) (- (+ (* 4 x1) (* 4 x2)) (* 8 x8)))
        A (+ (- (* (replace-div-by-zero 1 (* 192 (math/expt c2 2))) (math/expt (- (* 4 c1) c2) 4))
                (* (replace-div-by-zero 1 (* 24 c2)) (math/expt (- (* 4 c1) c2) 3))
                (* (replace-div-by-zero 1 (* 3840 (math/expt c2 3))) (math/expt (- (* 4 c1) c2) 5)))
             (* (replace-div-by-zero 1 (* 3840 (math/expt c2 3))) (math/expt (+ (* 4 c1) (* 3 c2)) 5)))
        B (+ (- (* (replace-div-by-zero 1 (* 192 (math/expt c2 2))) (math/expt (+ (* 4 c1) c2) 4)))
             (* (replace-div-by-zero 1 (* 24 c2)) (math/expt (+ (* 4 c1) c2) 3))
             (* (replace-div-by-zero 1 (* 3840 (math/expt c2 3))) (math/expt (+ (* 4 c1) c2) 5))
             (- (* (replace-div-by-zero 1 3840) (math/expt c2 3) (math/expt (- (* 4 c1) (* 3 c2)) 5))))
        C (* (replace-div-by-zero 4 15) (math/expt c2 2))
        D (+ (* (replace-div-by-zero 1 (* 192 (math/expt c2 2))) (math/expt (- (* 4 c2) c1) 4))
             (- (* (replace-div-by-zero 1 (* 3840 (math/expt c2 3))) (math/expt (- (* 4 c2) c1) 5)))
             (* (replace-div-by-zero 1 (* 7680 (math/expt c2 3))) (math/expt (+ (* 4 c2) (* 8 c1)) 5))
             (- (* (replace-div-by-zero 7 (* 7680 (math/expt c2 3))) (math/expt (- (* 4 c2) (* 8 c1)) 5)))
             (* (replace-div-by-zero 1 (* 768 (math/expt c2 3))) (math/expt (- (* 8 c1)) 5))
             (- (* (replace-div-by-zero c1 (* 96 (math/expt c2 3))) (math/expt (- (* 4 c2) (* 8 c1)) 4)))
             (* (replace-div-by-zero (- (* 2 c1) 1) (* 192 (math/expt c2 3))) (- (math/expt (* 8 c1) 4))))
        E (+ (* (/ 8 3) (math/expt c1 2)) (* (/ 1 30) (math/expt c2 2)))
        F (- (* (/ 2 3) c1 c2) (* (/ 1 30) (math/expt c2 2)))
        G (- (- (* (/ 16 3) (math/expt c1 2))) (* (/ 4 3) c1 c2) (* (/ 2 15) (math/expt c2 2)))
        H (+ (* (/ 2 3) c1 c2) (* (/ 1 30) (math/expt c2 2)))
        I (- (- (* (/ 16 3) (math/expt c1 2))) (* (/ 2 3) (math/expt c2 2)))
        J (* (/ 2 15) (math/expt c2 2))
        K (- (* (/ 4 3) c1 c2))
        u (m/matrix [[A E 0 0 (- F) 0 (- F) G F F]
                     [E B 0 0 (- H) 0 (- H) I H H]
                     [0 0 0 0 0 0 0 0 0 0]
                     [0 0 0 0 0 0 0 0 0 0]
                     [(- F) (- H) 0 0 C 0 J (- K) (- C) (- J)]
                     [0 0 0 0 0 0 0 0 0 0]
                     [(- F) (- H) 0 0 J 0 C (- K) (- J) (- C)]
                     [G I 0 0 (- K) 0 (- K) D K K]
                     [F H 0 0 (- C) 0 (- J) K C J]
                     [F H 0 0 (- J) 0 (- C) K J C]])
        temp (m/zero-matrix 30 30)]
    (->
     (->> temp
          (insert-submatrix-in-matrix u {:i 0 :j 0})
          (insert-submatrix-in-matrix u {:i 10 :j 10})
          (insert-submatrix-in-matrix u {:i 20 :j 20}))
     (m/mmul (* EI jacob)))))

(defn create-local-b
  "Create Local b system"
  [element f]
  (let [t (m/matrix [[59] [-1] [-1] [-1]
                     [4] [4] [4] [4] [4] [4]])
        J (jacobian element)
        Nt (m/zero-matrix 30 3)
        f-mat (m/matrix [[(-> f :x)]
                         [(-> f :y)]
                         [(-> f :z)]])]
    (->
     (->> Nt
          (insert-submatrix-in-matrix t {:i 0 :j 0})
          (insert-submatrix-in-matrix t {:i 10 :j 1})
          (insert-submatrix-in-matrix t {:i 20 :j 2})
          (m/mul (/ J 120.)))
     (m/mmul f-mat))))

(defn parse-params
  "Parses params from <buffer>"
  [buffer]
  (let [line (first (line-seq buffer))]
    (->> line
         (re-seq #"[^ ]+")
         (map #(Float/parseFloat %))
         (doall))))

(defn parse-quantities
  "Parse quantities from <buffer>"
  [buffer]
  (let [line (first (line-seq buffer))]
    (->> line
         (re-seq #"[^ ]+")
         (map #(Integer/parseInt %))
         (doall))))

(defn parse-nodes
  "Parse <quantity> nodes from <buffer>"
  [buffer nodes-quantity]
  (->> buffer
       (line-seq)
       (drop 2)
       (take nodes-quantity)
       (map #(re-seq #"[^ ]+" %))
       (map (fn [line] (map #(Float/parseFloat %) line)))
       (map (fn [[id x y z]] (->Node (int id) x y z)))
       (doall)))

(defn parse-elements
  "Parse <quantity> elements from <buffer>"
  [buffer elems-quantity nodes]
  (->> buffer
       (line-seq)
       (drop 3)
       (take elems-quantity)
       (map #(re-seq #"[^ ]+" %))
       (map (fn [line] (map #(Integer/parseInt %) line)))
       (map (fn [[id n1 n2 n3 n4 n5 n6 n7 n8 n9 n10]]
              (->Element id (get-node-by-id n1 nodes) (get-node-by-id n2 nodes)
                         (get-node-by-id n3 nodes) (get-node-by-id n4 nodes)
                         (get-node-by-id n5 nodes) (get-node-by-id n6 nodes)
                         (get-node-by-id n7 nodes) (get-node-by-id n8 nodes)
                         (get-node-by-id n9 nodes) (get-node-by-id n10 nodes))))
       (doall)))

(defn parse-condition-values
  "Get condition values"
  [buffer quantity nodes]
  (->> buffer
       (line-seq)
       (drop 3)
       (take quantity)
       (map #(re-seq #"[^ ]+" %))
       (map (fn [line] (map #(Float/parseFloat %) line)))
       (map (fn [[node-id value]]
              (->Condition
               (get-node-by-id (int node-id) nodes)
               value)))
       (doall)))

(defn read-input
  "Read from file"
  [file]
  (with-open [buffer (clojure.java.io/reader file)]
    (let [[EI fx fy fz] (parse-params buffer)
          [nodes-len elements-len dirichletx-len
           dirichlety-len dirichletz-len neumman-len] (parse-quantities buffer)
          nodes (parse-nodes buffer nodes-len)
          elements (parse-elements buffer elements-len nodes)]
      {:params {:EI EI, :f {:x fx :y fy :z fz}}
       :quantities {:nodes nodes-len
                    :elements elements-len
                    :dirichlet {:x dirichletx-len
                                :y dirichlety-len
                                :z dirichletz-len}
                    :neumann neumman-len}
       :nodes nodes
       :elements elements
       :conditions {:dirichlet {:x (parse-condition-values buffer dirichletx-len nodes)
                                :y (parse-condition-values buffer dirichlety-len nodes)
                                :z (parse-condition-values buffer dirichletz-len nodes)}
                    :neumann (parse-condition-values buffer neumman-len nodes)}})))

(defn assembly
  "returns the assemblied K and b from local systems"
  [nodes elements local-k-array local-b-array]
  (let [nodes-quant (count nodes)
        element-nodes (get-nodes-by-element (first elements))
        K (m/zero-matrix nodes-quant nodes-quant)
        b (m/zero-matrix nodes-quant 1)]
    (doseq [[index elem] (map-indexed vector elements)]
      (let [local-K (nth local-k-array index)
            local-b (nth local-b-array index)]
        (doseq [i (range 1 element-nodes)
                j (range 1 element-nodes)]
          (m/mset! K (node-index (get-in elem [(keyword (str "n" i))])) (node-index (get-in elem [(keyword (str "n" j))]))
                   (+ (m/mget K (node-index (get-in elem [(keyword (str "n" i))])) (node-index (get-in elem [(keyword (str "n" j))])))
                      (m/mget local-K (- i 1) (- j 1)))))
        (doseq [i (range 1 element-nodes)]
          (m/mset! b (node-index (get-in elem [(keyword (str "n" i))])) 0
                   (+ (m/mget b (node-index (get-in elem [(keyword (str "n" i))])) 0)
                      (m/mget local-b (- i 1) 0))))))
    {:K K
     :b b}))

(defn apply-neumann-condition
  "Applies neumann condition to b, returns b altered"
  [conditions b]
  (let [temp (m/matrix b)]
    (doseq [[condition] (map vector conditions)]
      (let [node (get-in condition [:node])
            value (get-in condition [:value])]
        (m/mset! temp (node-index node) 0
                 (+ (m/mget temp (node-index node) 0) value))))
    temp))

(defn dirichlet-row-helper
  "Removes rows"
  [K b indexes]
  (if (empty? indexes)
    {:K K
     :b b}
    (let [index (first indexes)
          new-K (remove-row K index)
          new-b (remove-row b index)]
      (dirichlet-row-helper new-K new-b (rest indexes)))))

(defn dirichlet-cols-helper
  "Removes cols"
  [K indexes]
  (if (empty? indexes)
    {:K K}
    (let [index (first indexes)
          new-K (remove-col K index)]
      (dirichlet-cols-helper new-K (rest indexes)))))

(defn apply-dirichlet-condition
  "Applies dirichlet condition to K, returns K altered"
  [conditions K b]
  (let [index-to-remove (map-indexed (fn [i condition]
                                       (- (node-index (get-in condition [:node])) i)) conditions)
        {temp-K :K new-b :b} (dirichlet-row-helper K b index-to-remove)
        {new-K :K} (dirichlet-cols-helper temp-K index-to-remove)]
    {:K new-K
     :b new-b}))

(defn calculate-fem
  "Last step"
  [K b]
  (let [K-cleand (m/emap #(if (zero? %)
                            (rand-int 2)
                            %) K)
        K-inv (m/inverse K-cleand)]
    (m/mmul K-inv b)))

(defn insert-row-helper
  "Insert dirichlet conditions to answers"
  [X conditions]
  (if (empty? conditions)
    X
    (let [condition (first conditions)]
      (insert-row-helper
       (add-row X (node-index (-> condition :node))
                [(-> condition :value)])
       (rest conditions)))))

(defn insert-condition-values
  "Inserts condition values in X vector"
  [X conditions]
  (insert-row-helper X conditions))

(defn write-to-file
  "writes X to file" [file X]
  (with-open [wrtr (clojure.java.io/writer file)]
    (.write wrtr (str "GiD Post Results File 1.0\n"))
    (.write wrtr (str "Result \"X\" \"Load Case 1\" 1 Scalar OnNodes\n"))
    (.write wrtr (str "ComponentNames \"X\"\n"))
    (.write wrtr (str "Values\n"))
    (doseq [[i x] (map-indexed vector X)]
      (.write wrtr (str (+ i 1) "\t" (m/mget x 0) "\n")))
    (.write wrtr (str "End Values"))))
