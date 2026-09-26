
  public:

    /// re-optimizes from the flow and the potentials of the last run
    /** Re-optimizes from the flow and the potentials of the last run instead
     * of starting from the zero flow and the zero potentials, which is what
     * a primal-dual algorithm allows after any change of the costs, of the
     * capacities or of the supplies: the flow of each arc is put within its
     * new bounds, and then at the bound that the reduced cost of the arc,
     * taken with the potentials of the last run, says (the upper one if the
     * reduced cost is negative, the lower one if it is positive), so that
     * the reduced cost of every residual arc is nonnegative again; what this
     * leaves at the nodes is an excess, which the shortest augmenting paths
     * of start() send where it is due, as they do in a run from scratch. It
     * does what run() does if the last run did not end at an optimal flow,
     * with lower bounds, if the supplies of the last run or the new ones do
     * not sum to zero, the root being set up differently then [see init()],
     * and if an arc of infinite capacity has a negative reduced cost, since
     * its flow cannot be put at its upper bound. */

    ProblemType runWarm(int factor = 4) {
      if (!_warm || _has_lower || _sum_supply != 0) return runCold(factor);
      Value sum = 0;
      for (int i = 0; i != _root; ++i) sum += _supply[i];
      if (sum != 0) return runCold(factor);

      const Value MAX = std::numeric_limits<Value>::max();
      _factor = factor;
      for (int i = 0; i != _root; ++i) _excess[i] = _supply[i];
      _excess[_root] = 0;

      // the flow of the last run, within the new bounds and at the one the
      // reduced cost says; the last arc out of a node is the one to the root
      for (int u = 0; u != _root; ++u) {
        const int last_out = _first_out[u + 1] - 1;
        for (int j = _first_out[u]; j != last_out; ++j) {
          if (!_forward[j]) continue;
          const int v = _target[j], rj = _reverse[j];
          const Value up = _upper[j];
          Value x = _res_cap[rj];
          if (x > up) x = up;
          const Cost c = _cost[j] + _pi[u] - _pi[v];
          if (c < 0) {
            if (up >= MAX) return runCold(factor);
            x = up;
          } else if (c > 0) {
            x = 0;
          }
          _res_cap[j] = up >= MAX ? INF : up - x;
          _res_cap[rj] = x;
          _excess[u] -= x;
          _excess[v] += x;
        }
      }
      for (int a = _first_out[_root]; a != _res_arc_num; ++a) {
        _res_cap[a] = 1;
        _res_cap[_reverse[a]] = 0;
      }

      // the scaling starts from the largest excess, as in init()
      if (_factor > 1) {
        Value max_sup = 0, max_dem = 0, max_cap = 0;
        for (int i = 0; i != _root; ++i) {
          const Value ex = _excess[i];
          if (ex > max_sup) max_sup = ex;
          if (-ex > max_dem) max_dem = -ex;
          const int last_out = _first_out[i + 1] - 1;
          for (int j = _first_out[i]; j != last_out; ++j)
            if (_res_cap[j] > max_cap) max_cap = _res_cap[j];
        }
        max_sup = std::min(std::min(max_sup, max_dem), max_cap);
        for (_delta = 1; 2 * _delta <= max_sup; _delta *= 2) ;
      } else {
        _delta = 1;
      }

      const ProblemType pt = start();
      _warm = (pt == OPTIMAL);
      return pt;
    }

  private:

    // run(), after which the flow and the potentials are there to reuse if
    // it has ended at an optimal flow [see runWarm()]
    ProblemType runCold(int factor) {
      const ProblemType pt = run(factor);
      _warm = (pt == OPTIMAL);
      return pt;
    }

    bool _warm = false;

