
  public:

    /// re-optimizes from the flow of the last run
    /** Re-optimizes from the flow of the last run instead of from a feasible
     * flow found anew: a cycle-canceling algorithm needs only a feasible
     * flow to start from, which the flow of the last run is after a change
     * of the costs, and which after a change of the capacities or of the
     * supplies it becomes once put within the new bounds and once what this
     * leaves at the nodes is sent where it is due [see feasibleFlow()], so
     * that the negative cycles are canceled from there. It does what run()
     * does if the last run did not end at an optimal flow, with lower
     * bounds, with an infinite upper bound, which init() bounds by the
     * deficits, when the supplies do not sum to zero and when no flow meets
     * the new bounds and supplies, which run() then proves. */

    ProblemType runWarm(Method method = CANCEL_AND_TIGHTEN) {
      if (!_warm || _has_lower || _sum_supply != 0) return runCold(method);

      const Value MAX = std::numeric_limits<Value>::max();
      Value sum = 0;
      for (int i = 0; i != _root; ++i) sum += _supply[i];
      if (sum != 0) return runCold(method);
      for (int j = 0; j != _res_arc_num; ++j)
        if (_forward[j] && _source[j] != _root && _upper[j] >= MAX)
          return runCold(method);
      int sup_num;
      if (!feasibleFlow(true, Value(1), sup_num)) return runCold(method);

      // the zero potentials init() would start from
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

