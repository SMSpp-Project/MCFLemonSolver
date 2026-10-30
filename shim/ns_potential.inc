    // Update potentials in the subtree that has been moved, or in the rest
    // of the tree if that is the smaller part: the rest follows the subtree
    // in the thread, up to the first node of the subtree, and its potentials
    // go the other way, which leaves every reduced cost as it is and moves
    // the potential of the root away from zero, by no more than _pi_lim,
    // the potentials being shifted back at the end of start()
    void updatePotential() {
      Cost sigma = _pi[v_in] - _pi[u_in] -
                   _pred_dir[u_in] * _cost[in_arc];
      int beg = u_in, end = _thread[_last_succ[u_in]];
      if (2 * _succ_num[u_in] > _node_num + 1) {
        const Cost r = _pi[_root];
        if (sigma >= r - _pi_lim && sigma <= r + _pi_lim) {
          beg = end;
          end = u_in;
          sigma = -sigma;
        }
      }
      for (int u = beg; u != end; u = _thread[u]) {
        _pi[u] += sigma;
      }
    }

    // how far the potential of the root may move from zero [see
    // ns_pi_lim()], set with the cost of the artificial arcs
    Cost _pi_lim = Cost(0);
