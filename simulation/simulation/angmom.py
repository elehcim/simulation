from pynbody.analysis.angmom import calc_faceon_matrix, calc_sideon_matrix


def sideon(h, vec_to_xform=calc_sideon_matrix,
             cen=None, vcen=None, move_all=True, **kwargs):
    """
    Reposition and rotate the simulation containing the halo h to see
    h's disk edge on.

    Given a simulation and a subview of that simulation (probably the
    halo of interest), this routine centers the simulation and rotates
    it so that the disk lies in the x-z plane. This gives a side-on
    view for SPH images, for instance.

    """

    from pynbody.analysis.angmom import ang_mom_vec
    from pynbody.analysis import halo
    from pynbody import filt, units, transformation

    if move_all:
        top = h.ancestor
    else:
        top = h

    # Top is the top-level view of the simulation, which will be
    # transformed

    if cen is None:
        # logger.info("Finding halo center...")
        # or h['pos'][h['phi'].argmin()]
        cen = halo.center(h, retcen=True, **kwargs)
        # logger.info("... cen=%s" % cen)

    tx = transformation.inverse_translate(top, cen)

    if vcen is None:
        vcen = halo.vel_center(h, retcen=True)

    tx = transformation.inverse_v_translate(tx, vcen)

    # logger.info("Calculating angular momentum vector...")
    am = ang_mom_vec(h)
    # print(" before rotation L:", am)
    trans = vec_to_xform(am)


    # logger.info("Transforming simulation...")

    tx = transformation.transform(tx, trans)
    # print(" after rotation  L:", ang_mom_vec(h))

    # logger.info("...done!")

    return tx


def faceon(h, **kwargs):
    return mysideon(h, vec_to_xform=calc_faceon_matrix, **kwargs)
