

function collapse_adjacent_leaves!(net::HybridNetwork)

    runagain = true
    while runagain
        runagain = false
        for leaf in net.leaf
            par = getparent(leaf)
            sib = getchildren(par)[1] == leaf ? getchildren(par)[2] : getchildren(par)[1]
            if sib.leaf && sib.name == leaf.name
                runagain = true
                PhyloNetworks.deleteleaf!(net, sib)
                break
            end
        end
    end

end

