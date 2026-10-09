
  public:

    /// re-optimizes from the flow of the last run
    /** Re-optimizes from the flow of the last run instead of from a feasible
     * flow found anew: a cycle-canceling algorithm needs only a feasible
     * flow to start from, and after a change of the costs the one of the
     * last run still is, so that init(), i.e., the circulation that finds a
     * feasible flow, is skipped and the negative cycles are canceled from
     * there. It does what run() does if the last run did not end at an
     * optimal flow, if that flow is no longer feasible (after a change of
     * the capacities or of the supplies), with lower bounds and when the
     * supplies do not sum to zero. */

    ProblemType runWarm(Method method = CANCEL_AND_TIGHTEN) {
      if (!_warm || _has_lower || _sum_supply != 0) return runCold(method);

      const Value MAX = std::numeric_limits<Value>::max();

      // the flow of the last run has to fit the bounds and the supplies
      ValueVector ex(_res_node_num, 0);
      for (int i = 0; i != _root; ++i) ex[i] = _supply[i];
      for (int j = 0; j != _res_arc_num; ++j) {
        if (!_forward[j] || _source[j] == _root || _target[j] == _root)
          continue;
        const Value x = _res_cap[_reverse[j]];
        if (_upper[j] >= MAX || x > _upper[j]) return runCold(method);
        ex[_source[j]] -= x;
        ex[_target[j]] += x;
      }
      for (int i = 0; i != _res_node_num; ++i)
        if (ex[i] != 0) return runCold(method);

      // the residual capacities under the new upper bounds, and the zero
      // potentials init() would start from
      for (int j = 0; j != _res_arc_num; ++j)
        if (_forward[j] && _source[j] != _root && _target[j] != _root)
          _res_cap[j] = _upper[j] - _res_cap[_reverse[j]];
      for (int i = 0; i != _res_node_num; ++i) _pi[i] = 0;

      start(method);
      _warm = true;
      return OPTIMAL;
    }

  private:

    // run(), after which the flow is there to reuse if it has ended at an
    // optimal one [see runWarm()]
    ProblemType runCold(Method method) {
      const ProblemType pt = run(method);
      _warm = (pt == OPTIMAL);
      return pt;
    }

    bool _warm = false;

