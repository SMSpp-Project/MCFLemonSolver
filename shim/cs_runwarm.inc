
  public:

    /// re-optimizes from the flow and the potentials of the last run
    /** Re-optimizes from the flow and the potentials of the last run instead
     * of starting from a flow found anew and the zero potentials: if the
     * flow of the last run is still feasible, i.e., after a change of the
     * costs, it is kept, the potentials are scaled back to the costs the
     * algorithm works on, and the scaling starts from the largest amount by
     * which a residual arc violates the optimality conditions under them,
     * rather than from the largest cost (the first phase being the one after
     * it, since the pair is already that close to optimal), so that a small
     * change of the costs costs a few phases at small epsilon. It does what run() does if the
     * last run did not end at an optimal flow, if the flow is no longer
     * feasible (after a change of the capacities or of the supplies, which
     * a feasible flow has to be found for anew [see init()]), with lower
     * bounds and when the supplies do not sum to zero. */

    ProblemType runWarm(Method method = PARTIAL_AUGMENT, int factor = 16) {
      if (!_warm || _has_lower || _sum_supply != 0)
        return runCold(method, factor);

      const Value MAX = std::numeric_limits<Value>::max();

      // the flow of the last run has to fit the bounds and the supplies
      ValueVector ex(_res_node_num, 0);
      for (int i = 0; i != _root; ++i) ex[i] = _supply[i];
      for (int j = 0; j != _res_arc_num; ++j) {
        if (!_forward[j] || _source[j] == _root || _target[j] == _root)
          continue;
        const Value x = _res_cap[_reverse[j]];
        if (_upper[j] >= MAX || x > _upper[j]) return runCold(method, factor);
        ex[_source[j]] -= x;
        ex[_target[j]] += x;
      }
      for (int i = 0; i != _res_node_num; ++i)
        if (ex[i] != 0) return runCold(method, factor);

      _alpha = factor;
      const LargeCost scale = static_cast<LargeCost>(_res_node_num) * _alpha;
      for (int i = 0; i != _res_node_num; ++i) {
        _excess[i] = 0;
        _pi[i] *= scale;
      }

      // the residual capacities under the new upper bounds, the costs the
      // algorithm works on, and how far the potentials are from optimal
      _epsilon = 0;
      for (int j = 0; j != _res_arc_num; ++j) {
        const int u = _source[j], v = _target[j];
        if (u == _root || v == _root) continue;
        if (_forward[j]) _res_cap[j] = _upper[j] - _res_cap[_reverse[j]];
        _cost[j] = static_cast<LargeCost>(_scost[j]) * scale;
      }
      for (int j = 0; j != _res_arc_num; ++j) {
        if (_res_cap[j] <= 0) continue;
        const LargeCost rc = _cost[j] + _pi[_source[j]] - _pi[_target[j]];
        if (-rc > _epsilon) _epsilon = -rc;
      }

      // the pair is already epsilon-optimal, hence the first phase is the
      // one after epsilon, as the loop of start() would take it, rather than
      // a refinement to where the pair already is
      _epsilon = _epsilon >= _alpha ? _epsilon / _alpha
                                    : (_epsilon >= 1 ? 1 : _epsilon);

      _max_rank = _alpha * _res_node_num;
      _buckets.resize(_max_rank);

      start(method);
      _warm = true;
      return OPTIMAL;
    }

  private:

    // run(), after which the flow and the potentials are there to reuse if
    // it has ended at an optimal flow [see runWarm()]
    ProblemType runCold(Method method, int factor) {
      const ProblemType pt = run(method, factor);
      _warm = (pt == OPTIMAL);
      return pt;
    }

    bool _warm = false;

