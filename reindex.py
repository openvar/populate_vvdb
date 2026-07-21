from VariantValidator import update_vv_db

# Re index
db = update_vv_db.connect()
update_vv_db.drop_core_indexes(db)
update_vv_db.rebuild_core_indexes(db)