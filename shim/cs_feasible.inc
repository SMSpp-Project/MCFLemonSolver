
    // finds a flow within the bounds that meets the supplies, which init()
    // of LEMON finds with Circulation on the graph and its maps: here with
    // the same push-relabel algorithm, i.e., the active node of highest
    // level first, every node at level zero at the start, a push towards
    // any lower level and the first empty level as the proof that no such
    // flow exists, but on the residual network the algorithm works on, and
    // from either the flow the arcs into the nodes with a deficit take
    // greedily, as Circulation does, or, with warm, the flow of the last run
    // put within the new bounds, so that after a change of the capacities
    // or of the supplies only what the change has displaced is sent anew
    // [see runWarm()]. It sets the residual capacities of the arcs of the
    // graph and of those of the root, the ones out of the root having
    // root_cap if the supplies sum to zero, as init() does, and gives in
    // sup_num the number of nodes whose supply net of the lower bounds is
    // positive; it returns false if there is no such flow.

    bool feasibleFlow(bool warm, Value root_cap, int& sup_num) {
      const int n = _root;
      ValueVector ex(_supply.begin(), _supply.begin() + n);
      if (_has_lower) {
        for (int u = 0; u != n; ++u) {
          const int last_out = _first_out[u + 1] - 1;
          for (int j = _first_out[u]; j != last_out; ++j) {
            if (!_forward[j]) continue;
            ex[u] -= _lower[j];
            ex[_target[j]] += _lower[j];
          }
        }
      }
      sup_num = 0;
      for (int u = 0; u != n; ++u)
        if (ex[u] > 0) ++sup_num;

      // the flow to start from; the last arc out of a node is the one to
      // the root, which is not one of the graph
      for (int u = 0; u != n; ++u) {
        const int last_out = _first_out[u + 1] - 1;
        for (int j = _first_out[u]; j != last_out; ++j) {
          if (!_forward[j]) continue;
          const int v = _target[j], rj = _reverse[j];
          const Value cap = _has_lower ? _upper[j] - _lower[j] : _upper[j];
          Value f = warm ? _res_cap[rj] : -ex[v];
          if (f > cap) f = cap;
          if (f < 0) f = 0;
          _res_cap[j] = cap - f;
          _res_cap[rj] = f;
          ex[u] -= f;
          ex[v] += f;
        }
      }

      // the active nodes, i.e., those with a positive excess, in a list for
      // each level, and how many nodes each level has
      IntVector level(n, 0), count(n + 1, 0), next(n, -1), head(n + 1, -1);
      count[0] = n;
      for (int u = 0; u != n; ++u)
        if (ex[u] > 0) {
          next[u] = head[0];
          head[0] = u;
        }
      for (int top = 0; ; ) {
        while (top >= 0 && head[top] < 0) --top;
        if (top < 0) break;
        const int u = head[top];
        Value e = ex[u];
        int low = n;
        const int last_out = _first_out[u + 1] - 1;
        for (int j = _first_out[u]; j != last_out; ++j) {
          const Value r = _res_cap[j];
          const int v = _target[j];
          if (r <= 0 || v == u) continue;
          const int lv = level[v];
          if (lv >= top) {
            if (lv < low) low = lv;
            continue;
          }
          const Value d = r < e ? r : e;
          _res_cap[j] = r - d;
          _res_cap[_reverse[j]] += d;
          if (ex[v] <= 0 && ex[v] + d > 0) {
            next[v] = head[lv];
            head[lv] = v;
          }
          ex[v] += d;
          if ((e -= d) == 0) break;
        }
        ex[u] = e;
        head[top] = next[u];
        if (e == 0) continue;
        // u is lifted above the lowest level it has a residual arc to, and
        // no such arc, or no node left on its level, proves that no flow
        // exists [see the barrier of Circulation]
        if (low == n || --count[top] == 0) return false;
        level[u] = ++low;
        ++count[low];
        next[u] = head[low];
        head[low] = u;
        top = low;
      }

      // the arcs of the root: with supplies that sum to less than zero, the
      // ones out of it carry what is left of the deficits, as in init()
      for (int a = _first_out[_root]; a != _res_arc_num; ++a) {
        const int ra = _reverse[a];
        if (_sum_supply < 0) {
          _res_cap[a] = -_sum_supply + 1;
          _res_cap[ra] = -ex[_target[a]];
        } else {
          _res_cap[a] = root_cap;
          _res_cap[ra] = 0;
        }
        _cost[a] = 0;
        _cost[ra] = 0;
      }
      return true;
    }
