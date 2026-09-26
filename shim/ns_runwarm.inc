      // the tail of run(): a basis is there to reuse if this run has ended
      // at an optimal one, or at one that proves the instance infeasible,
      // which is optimal for the problem with the artificial arcs and makes
      // a start after a change as good as any [see runWarm()]; not with
      // lower bounds, which start() does not put back in the supplies when
      // it answers INFEASIBLE
      _repair = false;
      ProblemType pt = start(pivot_rule);
      _warm = (pt == OPTIMAL || (pt == INFEASIBLE && !_has_lower));
      return pt;
      }

    /// re-optimizes from the basis of the last run
    /** Re-optimizes from the basis of the last run instead of starting from
     * the artificial one. After a change of the costs alone the basis is
     * still primal feasible, hence the potentials are recomputed on its tree
     * and the simplex goes on from there, which is what a Lagrangian or a
     * Frank-Wolfe decomposition asks for. After a change of the capacities
     * or of the supplies the flow of the tree is recomputed first, and the
     * subtree below an arc whose flow falls outside its bounds is hung from
     * the root by its artificial arc [see repairBasis()], as init() hangs
     * every node, so that the basis is again primal feasible and the
     * artificial arcs are driven out by their cost. A change of the graph,
     * of the lower bounds or of the supplies that leaves them not summing to
     * zero drops the basis, and then this method does what run() does. */

    ProblemType runWarm(PivotRule pivot_rule = BLOCK_SEARCH) {
      if (!_warm) return run(pivot_rule);
      if (_repair) {
        if (!repairBasis()) return run(pivot_rule);
        _repair = false;
      }

      // the artificial arcs keep a cost larger than any of the real ones,
      // which the change may have made larger than it was
      Cost art = 0;
      if (std::numeric_limits<Cost>::is_exact) {
        art = std::numeric_limits<Cost>::max() / 2 + 1;
      } else {
        for (int i = 0; i != _arc_num; ++i)
          if (_cost[i] > art) art = _cost[i];
        art = (art + 1) * _node_num;
      }
      for (int e = _arc_num; e != _all_arc_num; ++e)
        if (_cost[e] != 0) _cost[e] = art;

      // the potentials of the basis tree, in the order of the thread, which
      // gives every node after its own father
      _pi[_root] = 0;
      for (int u = _thread[_root]; u != _root; u = _thread[u]) {
        int e = _pred[u];
        _pi[u] = _pred_dir[u] == DIR_UP ? _pi[_parent[u]] - _cost[e]
                                        : _pi[_parent[u]] + _cost[e];
      }

      ProblemType ptw = start(pivot_rule);
      _warm = (ptw == OPTIMAL || (ptw == INFEASIBLE && !_has_lower));
      return ptw;
      }

    /// takes in the new upper bound of an arc, keeping the basis if it can
    /** Takes in the new upper bound of the arc of internal index i, as
     * upperMap() gives it. The basis of the last run is kept as it is if the
     * flow the arc holds stays inside the new bound and the arc does not sit
     * at its upper bound, which is the case of the arc that is closed while
     * it carries no flow and of the one that is opened again; otherwise the
     * flow of the tree has to be recomputed, which runWarm() does before
     * starting [see repairBasis()]. With lower bounds the basis is dropped
     * instead, the capacity the algorithm works on being that of the problem
     * the lower bounds have been taken out of, as init() computes it. */

    void warmUpper(int i, Value up) {
      if (up == _upper[i]) return;
      _upper[i] = up;
      if (!_warm) return;

      Value cap;
      if (_has_lower) {
        Value c = _lower[i];
        if (c >= 0) {
          cap = up < MAX ? up - c : INF;
        } else {
          cap = up < MAX + c ? up - c : INF;
        }
      } else {
        cap = up;
      }

      const bool moves = _state[i] == STATE_UPPER || _flow[i] > cap;
      if (moves && _has_lower) { _warm = false; return; }
      _cap[i] = cap;
      if (moves) _repair = true;
      }

    /// takes in the new supply of a node, keeping the basis
    /** Takes in the new supply of the node of internal index i, as
     * supplyMap() gives it: the basis of the last run is kept, and the flow
     * of its tree is recomputed by runWarm() before starting [see
     * repairBasis()]. */

    void warmSupply(int i, Value s) {
      if (s == _supply[i]) return;
      _supply[i] = s;
      if (_warm) _repair = true;
      }

    /// makes the basis of the last run primal feasible again
    /** Makes the basis of the last run primal feasible again after a change
     * of the capacities or of the supplies, keeping its tree: every arc out
     * of the tree is put at the bound it sits at, and the flow of the tree
     * arcs is recomputed from the leaves up, each carrying what the subtree
     * below it has to send. An arc whose flow falls outside its bounds leaves
     * the tree at the bound it violates, and the subtree below it is hung
     * from the root by the artificial arc of its top node, directed as the
     * flow it has to carry and with the cost init() would give it, i.e.,
     * zero towards the root and larger than any real one from it, so that
     * the simplex drives it out as it does from the artificial basis. The
     * thread and the sizes of the subtrees are then rebuilt from the fathers
     * by a depth-first visit. It returns false, leaving the basis to be
     * dropped, if there are lower bounds or if the supplies of the last run
     * or the new ones do not sum to zero, the artificial arcs being built
     * differently then [see init()]. */

    bool repairBasis() {
      if (_has_lower || _sum_supply != 0) return false;
      Value sum = 0;
      for (int u = 0; u != _node_num; ++u) sum += _supply[u];
      if (sum != 0) return false;

      // what each node has to send along the tree: its supply, net of the
      // flow of the arcs out of the tree, each at the bound it sits at
      std::vector<Value> net(_node_num + 1, 0);
      for (int u = 0; u != _node_num; ++u) net[u] = _supply[u];
      for (int e = 0; e != _arc_num; ++e) {
        if (_state[e] == STATE_TREE) continue;
        if (_state[e] == STATE_UPPER && _cap[e] >= MAX)
          _state[e] = STATE_LOWER;  // an upper bound no longer there
        _flow[e] = _state[e] == STATE_UPPER ? _cap[e] : 0;
        net[_source[e]] -= _flow[e];
        net[_target[e]] += _flow[e];
      }
      for (int e = _arc_num; e != _all_arc_num; ++e)
        if (_state[e] != STATE_TREE) _flow[e] = 0;

      // the flow of the tree arcs, from the leaves up, i.e., in the reverse
      // order of the thread, which gives every node after its subtree
      for (int u = _rev_thread[_root]; u != _root; u = _rev_thread[u]) {
        int e = _pred[u];
        Value n = net[u];
        if (e < _arc_num) {
          const int p = _parent[u];
          const Value f = _pred_dir[u] == DIR_UP ? n : -n;
          if (f >= 0 && f <= _cap[e]) {
            _flow[e] = f;
            net[p] += n;
            continue;
          }
          // the arc leaves the tree at the bound it violates, and what it
          // carries there goes from u to p
          _state[e] = f < 0 ? STATE_LOWER : STATE_UPPER;
          _flow[e] = f < 0 ? 0 : _cap[e];
          const Value s = _pred_dir[u] == DIR_UP ? _flow[e] : -_flow[e];
          n -= s;
          net[p] += s;
          e = _arc_num + u;
          _state[e] = STATE_TREE;
          _parent[u] = _root;
          _pred[u] = e;
        }
        // u hangs from the root by its artificial arc, directed as n says
        if (n >= 0) {
          _pred_dir[u] = DIR_UP;
          _source[e] = u;
          _target[e] = _root;
          _flow[e] = n;
          _cost[e] = 0;
        } else {
          _pred_dir[u] = DIR_DOWN;
          _source[e] = _root;
          _target[e] = u;
          _flow[e] = -n;
          _cost[e] = 1;  // nonzero, runWarm() gives it the artificial cost
        }
        net[_root] += n;
      }

      // the thread and its reverse, the sizes of the subtrees and their last
      // nodes, from the fathers, by a depth-first visit from the root
      std::vector<int> child(_node_num + 1, -1), sibling(_node_num + 1, -1);
      for (int u = _node_num - 1; u >= 0; --u) {
        sibling[u] = child[_parent[u]];
        child[_parent[u]] = u;
      }
      std::vector<int> order, stack(1, _root);
      order.reserve(_node_num + 1);
      while (!stack.empty()) {
        const int u = stack.back();
        stack.pop_back();
        order.push_back(u);
        for (int v = child[u]; v != -1; v = sibling[v]) stack.push_back(v);
      }
      const int k = static_cast<int>(order.size());
      if (k != _node_num + 1) return false;  // not a tree: cannot happen

      for (int i = 0; i != k; ++i) {
        const int u = order[i], v = order[(i + 1) % k];
        _thread[u] = v;
        _rev_thread[v] = u;
        _succ_num[u] = 1;
      }
      for (int i = k - 1; i > 0; --i)
        _succ_num[_parent[order[i]]] += _succ_num[order[i]];
      for (int i = 0; i != k; ++i)
        _last_succ[order[i]] = order[i + _succ_num[order[i]] - 1];

      return true;
      }

    /// the cycle that proves the instance unbounded
    /** The directed cycle of negative cost and infinite capacity that proves
     * the instance unbounded, which the algorithm has in its hands when it
     * says so: the cycle of the pivot it cannot perform, i.e., the entering
     * arc plus the path of the basis tree between its two ends, every arc of
     * which is traversed forward, an arc traversed backward having to hold an
     * infinite flow for the step to be infinite. The cycle is given in the
     * direction the flow takes, and it is empty if the last run did not end
     * with UNBOUNDED or if the cycle goes through an artificial arc, which is
     * no cycle of the instance. */

    std::vector<Arc> unbCycle() const {
      std::vector<Arc> cycle;
      if (_unb_arc < 0) return cycle;

      std::vector<int> ids;
      ids.push_back(_unb_arc);
      for (int u = v_in; u != join; u = _parent[u]) ids.push_back(_pred[u]);
      for (int u = u_in; u != join; u = _parent[u]) ids.push_back(_pred[u]);

      // an artificial arc is none of the instance
      for (std::size_t i = 0; i != ids.size(); ++i)
        if (ids[i] >= _arc_num) return std::vector<Arc>();

      // the internal numbering of the arcs is that of _arc_id, whose inverse
      // is taken here, once, over the arcs of the graph
      std::vector<Arc> arc_of(_arc_num, INVALID);
      for (ArcIt a(_graph); a != INVALID; ++a) arc_of[_arc_id[a]] = a;

      for (std::size_t i = 0; i != ids.size(); ++i) {
        if (arc_of[ids[i]] == INVALID) return std::vector<Arc>();
        cycle.push_back(arc_of[ids[i]]);
      }
      return cycle;
