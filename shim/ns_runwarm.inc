      // the tail of run(): a basis is there to reuse only if this run has
      // ended at an optimal one [see runWarm()]
      ProblemType pt = start(pivot_rule);
      _warm = (pt == OPTIMAL);
      return pt;
      }

    /// re-optimizes after a change of the costs alone
    /** Re-optimizes after a change of the costs alone: the basis of the last
     * run is still primal feasible, hence the potentials are recomputed on
     * its tree and the simplex goes on from there instead of starting from
     * the artificial basis, which is what a Lagrangian or a Frank-Wolfe
     * decomposition asks for, changing nothing but the costs at each
     * iteration. A change of the capacities, of the supplies or of the graph
     * drops the basis, since it is no longer feasible or no longer there at
     * all, and then this method does what run() does. */

    ProblemType runWarm(PivotRule pivot_rule = BLOCK_SEARCH) {
      if (!_warm) return run(pivot_rule);

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
      _warm = (ptw == OPTIMAL);
      return ptw;
      }

    /// takes in the new upper bound of an arc, keeping the basis if it can
    /** Takes in the new upper bound of the arc of internal index i, as
     * upperMap() gives it, and tells whether the basis of the last run
     * survives it: it does if the flow the arc holds stays inside the new
     * bound and, when the arc is out of the basis tree, if it sits at its
     * lower bound, where it stays whatever the upper bound becomes; an arc
     * that sits at its upper bound would have to move with it, and the flow
     * of the tree would have to be recomputed, which is what init() does.
     * This is the case of the arc that is closed while it carries no flow,
     * and of the one that is opened again, both of which a decomposition
     * that closes and opens arcs does at each iteration [see runWarm()].
     * The capacity the algorithm works on is that of the problem the lower
     * bounds have been taken out of, as init() computes it. */

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

      if (_state[i] == STATE_UPPER || _flow[i] > cap) _warm = false;
      else _cap[i] = cap;
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
